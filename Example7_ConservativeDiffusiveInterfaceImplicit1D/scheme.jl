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
#   "extrap"    pe*(1 - pe), pe = 2*phio - phioo, one solve per step (BDF2
#               with an extrapolated explicit term, the usual IMEX choice);
#               on the first step phioo = phio, so pe = phio
# Returns (x, phi(T), dt, nsteps, snapshots, snap_times, ok), x in order.
#==============================================================================#

using Finch
using DelimitedFiles

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
    ctx = init_finch("cdi1dimplicit"; verbose = verbose, solver = :bicgstab, solver_tol = 1e-12,
                     target = JuliaCPUTarget(out_dir = proj_dir, func_name = "cdi1dimplicit"))
    domain(1)
    mesh(MESH_LINE; elsperdim = nel)
    function_space(order = 1)

    phi   = variable("phi")     # unknown: the current Picard iterate
    phio  = variable("phio")    # previous time level
    phioo = variable("phioo")   # two time levels back
    phis  = variable("phis")    # lagged field for the sharpening term
    gx    = variable("gx")      # projected gradient of phis, for the nodal normal
    a0    = variable("a0")      # BDF coefficients over dt, as fields: BDF1 on
    a1    = variable("a1")      # the first step, switched to BDF2 by an update
    a2    = variable("a2")
    test_symbol("v")

    # The profile spans the domain, so the ends carry the two bulk values
    boundary(phi, 1, BC_DIRICHLET, 0.0)
    boundary(phi, 2, BC_DIRICHLET, 1.0)

    # Initial condition: the equilibrium profile centred at x0, in every level
    ic = (x, y, z, t) -> 0.5*(1 + tanh((x - X0)/(2*eps)))
    for f in (phi, phio, phioo, phis)
        initial(f, ic)
    end
    initial(a0, 1.0/dt); initial(a1, -1.0/dt); initial(a2, 0.0)

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

    # The time loop. Every entry is quasi-static; the stepper only sets dt.
    time_interval(T)
    time_stepper(STEP_EULER_IMPLICIT; dt = dt, store_trajectory = snapshots > 0,
                 snapshot_every = max(1, nsteps ÷ max(snapshots, 1)))

    # Weak form of one Picard iteration (residual = 0). Terms with phi go to
    # the matrix, the rest to the right-hand side:
    #   theta1*(a0*phi + a1*phio + a2*phioo, v)          time derivative (BDF)
    #   - theta2/2*(u*phi, grad v)                        advection, conservative half
    #   + theta2/2*(u*grad phi, v)                        advection, other half
    #   + theta3*(gamma*eps*grad phi, grad v)             diffusion
    #   - theta4*(gamma*S(phi, phis)*nhat(phis), grad v)  sharpening, linearized
    # The boundary terms from integrating by parts drop out because v vanishes
    # at the two Dirichlet ends.
    weak_form(gx, "gx*v - grad(phis)*v")
    weak_form(phi, "theta1*(a0*phi + a1*phio + a2*phioo)*v" *
                   " - 0.5*theta2*u*phi*grad(v)" *
                   " + 0.5*theta2*u*grad(phi)*v" *
                   " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
                   " - theta4*gamma*" * SHARPENING[form] * "*gx*pow(gx*gx + delta, -0.5)*grad(v)")

    first_lag = form == "extrap" ? "2*phio - phioo" : "phio"
    entries = Any[update(phis => first_lag), gx, phi]
    for it = 2:iters
        append!(entries, [update(phis => "phi"), gx, phi])
    end
    append!(entries, [update(phioo => "phio"), update(phio => "phi"),
                      update(a0 => "$(1.5/dt)"), update(a1 => "$(-2.0/dt)"), update(a2 => "$(0.5/dt)")])
    step_schedule(entries...)

    solve()
    ok = true
    try
        Finch.run_emitted_julia!(ctx.target, phi)
    catch err
        ok = false
        println("run failed: ", sprint(showerror, err))
    end

    x = ctx.grid.nodes[1, :]
    o = sortperm(x)
    snaps = Vector{Vector{Float64}}()
    times = Float64[]
    if ok && snapshots > 0
        results = joinpath(proj_dir, "results")
        _, metas = Finch.read_results_manifest(results)
        snaps = [Finch.load_snapshot_values(results, s)["phi"][o] for s in metas]
        times = [Float64(s["t"]) for s in metas]
    end
    return (x = x[o], phi = phi.values[1, o], dt = dt, nsteps = nsteps,
            snapshots = snaps, snap_times = times, ok = ok)
end
