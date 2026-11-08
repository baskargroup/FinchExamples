#==============================================================================#
# The semi-implicit 1D CDI scheme of main.jl as a function, so that main.jl,
# convergence.jl and sharpening_sweep.jl run the same code. Finch V2 builds a
# fresh context (and generated project) for each call.
#
#   run_case(CFL, iters; form, proj_dir, snapshots)
#
# iters Picard iterations per step: iteration 1 lags with phio, each later
# one with the previous iterate. form selects how the sharpening term
# gamma*(phi - phi^2)*nhat is split between the unknown phi and the lagged
# field phis (nhat always from phis):
#   "derived"   phi*(1 - phis)                 as derived (default)
#   "lagged"    phis*(1 - phis)                fully lagged, a right-hand side
#   "reversed"  phis*(1 - phi)                 the other factor implicit
#   "newton"    (1 - 2*phis)*phi + phis^2      Newton about phis
#   "extrap"    pe*(1 - pe), pe = 2*phio - phi(n-1), one solve per step (BDF2
#               with an extrapolated explicit term, the usual IMEX choice);
#               on the first step phi(n-1) = phio, so pe = phio
# Returns (x, phi(T), dt, nsteps, snapshots, snap_times, ok), x in order.
#
# Finch V2 (branch dev): the level two steps back is Finch's own history
# `old(phi, 2)` and the BDF coefficients are coefficients of t (BDF1 on the
# first step, BDF2 after), so only phio and phis remain as fields. The target
# (CPU Julia, or MLIR on the CPU/GPU) is chosen by FINCH_TARGET, see ../common.jl.
#==============================================================================#

include(joinpath(@__DIR__, "..", "common.jl"))   # Finch, FINCH_TARGET switch, helpers

# Physical parameters (same case as Example 5)
const NEL  = 200
const DX   = 1.0/NEL
const UVEL = 1.0                # advection velocity
const GAM  = 1.0                # sharpening parameter, usually max |u|
const EPS  = DX                 # interface thickness; eps/dx >= (u/gamma + 1)/2
                                # is what keeps phi within [0,1]
const X0   = 0.25               # initial interface position
const T    = 0.25

const SHARPENING = Dict(
    "derived" => "phi*(1 - phis)",
    "lagged"  => "phis*(1 - phis)",
    "newton"  => "((1 - 2*phis)*phi + phis*phis)",
    "reversed" => "phis*(1 - phi)",
    "extrap"  => "phis*(1 - phis)")

function run_case(CFL, iters; form = "derived", nel = NEL,
                  proj_dir = mktempdir(), snapshots = 20, verbose = false)
    dx = 1.0/nel
    eps = EPS*NEL/nel                   # eps = dx on any mesh
    nsteps = ceil(Int, T/(CFL*dx/UVEL) - 1e-9)   # dt rounded down so that T is reached exactly
    dt = T/nsteps
    mkpath(proj_dir)

    # Advection makes the matrix non-symmetric: BiCGSTAB
    ctx = init_finch("cdi1dimplicit"; verbose = verbose, finch_solver_kwargs(:bicgstab; tol = 1e-12)...,
                     target = finch_target("cdi1dimplicit", proj_dir))
    domain(1)
    mesh(MESH_LINE; elsperdim = nel)
    function_space(order = 1)

    phi   = variable("phi")     # unknown: the current Picard iterate; old(phi, 2) is two levels back
    phio  = variable("phio")    # previous time level (old(phi) would be the Picard iterate)
    phis  = variable("phis")    # lagged field for the sharpening term
    gx    = variable("gx")      # projected gradient of phis, for the nodal normal
    test_symbol("v")

    # The profile spans the domain, so the ends carry the two bulk values
    boundary(phi, 1, BC_DIRICHLET, 0.0)
    boundary(phi, 2, BC_DIRICHLET, 1.0)

    # Initial condition: the equilibrium profile centred at x0, in every level
    ic = (x, y, z, t) -> 0.5*(1 + tanh((x - X0)/(2*eps)))
    for f in (phi, phio, phis)
        initial(f, ic)
    end
    # BDF coefficients over dt: BDF1 on the first step (the solves of step n see
    # t = n*dt), BDF2 after
    coefficient("a0", (x, y, z, t) -> ifelse(t < 1.5*dt,  1.0/dt,  1.5/dt))
    coefficient("a1", (x, y, z, t) -> ifelse(t < 1.5*dt, -1.0/dt, -2.0/dt))
    coefficient("a2", (x, y, z, t) -> ifelse(t < 1.5*dt,  0.0,     0.5/dt))

    # Coefficients. The theta are the model coefficients of the derivation,
    # all one in dimensional form.
    coefficient("theta1", 1.0)
    coefficient("theta2", 1.0)
    coefficient("theta3", 1.0)
    coefficient("theta4", 1.0)
    coefficient("u", UVEL)
    coefficient("gamma", GAM)
    coefficient("epsilon", eps)
    coefficient("delta", 1e-12)     # keeps nhat finite where grad(phis) = 0

    # The time loop. Every entry is quasi-static; the stepper only sets the
    # step and keeps the history old(phi, 2).
    time_interval(T)
    time_stepper(STEP_BDF2; Nsteps = nsteps, store_trajectory = snapshots > 0,
                 snapshot_every = max(1, nsteps ÷ max(snapshots, 1)))

    # Weak form of one Picard iteration (residual = 0). Terms with phi go to
    # the matrix, the rest to the right-hand side:
    #   theta1*(a0*phi + a1*phio + a2*old(phi, 2), v)    time derivative (BDF)
    #   - theta2/2*(u*phi, grad v)                        advection, conservative half
    #   + theta2/2*(u*grad phi, v)                        advection, other half
    #   + theta3*(gamma*eps*grad phi, grad v)             diffusion
    #   - theta4*(gamma*S(phi, phis)*nhat(phis), grad v)  sharpening, linearized
    # The boundary terms from integrating by parts drop out because v vanishes
    # at the two Dirichlet ends.
    weak_form(gx, "gx*v - grad(phis)*v")
    weak_form(phi, "theta1*(a0*phi + a1*phio + a2*old(phi, 2))*v" *
                   " - 0.5*theta2*u*phi*grad(v)" *
                   " + 0.5*theta2*u*grad(phi)*v" *
                   " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
                   " - theta4*gamma*" * SHARPENING[form] * "*" * NDOTGV_NODAL_1D * "*grad(v)")

    first_lag = form == "extrap" ? "2*phio - old(phi, 2)" : "phio"
    entries = Any[update(phis => first_lag), gx, phi]
    for it = 2:iters
        append!(entries, [update(phis => "phi"), gx, phi])
    end
    push!(entries, update(phio => "phi"))
    step_schedule(entries...)

    solve()
    ok = true
    try
        finch_run!(ctx, phi)
    catch err
        ok = false
        println("run failed: ", sprint(showerror, err))
    end

    x = ctx.grid.nodes[1, :]
    o = sortperm(x)
    snaps = Vector{Vector{Float64}}()
    times = Float64[]
    phiT  = fill(NaN, length(x))
    if ok
        # the final state is the last entry of results/ (both targets write it)
        times, fields = read_snapshots(proj_dir)
        check_final_time(times, T)
        phiT = final_values(fields, "phi")[o]
        snaps = snapshots > 0 ? [fields["phi"][o, k] for k in eachindex(times)] : snaps
        snapshots > 0 || (times = Float64[])
    end
    return (x = x[o], phi = phiT, dt = dt, nsteps = nsteps,
            snapshots = snaps, snap_times = times, ok = ok)
end
