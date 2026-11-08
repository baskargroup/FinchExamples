#==============================================================================#
# The semi-implicit 2D CDI scheme of main.jl as a function, so that main.jl
# and sharpening_sweep.jl run the same code. Finch V2 builds a fresh context
# (and generated project) for each call.
#
#   run_case(CFL, iters; form, nel, proj_dir, snapshots)
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
# Returns (x, y, phi(T), dt, nsteps, snapshots, snap_times, mass, ok); mass
# is int phi after every step (index 1 is t = 0).
#==============================================================================#

using Finch
using DelimitedFiles

# Physical parameters (same case as Example 6)
const NEL  = 300
const UX   = 1.0                # advection velocity components
const UY   = 0.0
const UMAX = sqrt(UX^2 + UY^2)
const GAM  = UMAX               # sharpening parameter gamma (= max |u|)
const R    = 0.15               # circle radius
const XC0  = 0.3                # initial circle center
const YC0  = 0.5
const T    = 0.4

const SHARPENING = Dict(
    "derived"  => "phi*(1 - phis)",
    "lagged"   => "phis*(1 - phis)",
    "reversed" => "phis*(1 - phi)",
    "newton"   => "((1 - 2*phis)*phi + phis*phis)",
    "extrap"   => "phis*(1 - phis)")

function run_case(CFL, iters; form = "derived", nel = NEL,
                  proj_dir = mktempdir(), snapshots = 60, verbose = false)
    dx  = 1.0/nel
    eps = dx                            # interface thickness (eps/dx >= (u/gamma + 1)/2)
    nsteps = ceil(Int, T/(CFL*dx/UMAX) - 1e-9)   # dt rounded down so that T is reached exactly
    dt = T/nsteps
    mkpath(proj_dir)

    # The matrix changes every Picard iteration and is not symmetric
    # (advection): BiCGSTAB. A loose tolerance would show up directly as a
    # mass drift of that size.
    ctx = init_finch("cdi2dimplicit"; verbose = verbose, solver = :bicgstab, solver_tol = 1e-12,
                     target = JuliaCPUTarget(out_dir = proj_dir, func_name = "cdi2dimplicit"))
    domain(2)
    mesh(MESH_QUAD; elsperdim = nel)   # unit square
    function_space(order = 1)

    # Boundary ids: 1 = x==0 (left, inflow), 2 = x==1 (right, outflow),
    # 3 = y==0 and 4 = y==1 (bottom and top, u.n = 0)

    phi   = variable("phi")     # unknown: the current Picard iterate
    phio  = variable("phio")    # previous time level
    phioo = variable("phioo")   # two time levels back
    phis  = variable("phis")    # lagged field for the sharpening term
    gx    = variable("gx")      # projected gradient of phis, x component
    gy    = variable("gy")      # and y component
    a0    = variable("a0")      # BDF coefficients over dt, as fields: BDF1 on
    a1    = variable("a1")      # the first step, switched to BDF2 by an update
    a2    = variable("a2")
    test_symbol("v")

    # Inflow: the fluid entering at x = 0 carries phi = 0, so that value is imposed.
    boundary(phi, 1, BC_DIRICHLET, 0.0)

    # Initial condition: equilibrium tanh profile around the circle, in every level
    ic = (x, y, z, t) -> 0.5*(1 - tanh((sqrt((x - XC0)^2 + (y - YC0)^2) - R)/(2*eps)))
    for f in (phi, phio, phioo, phis)
        initial(f, ic)
    end
    initial(a0, 1.0/dt); initial(a1, -1.0/dt); initial(a2, 0.0)

    # Coefficients (theta as in the derivation, dimensional form)
    coefficient("theta1", 1.0)
    coefficient("theta2", 1.0)
    coefficient("theta3", 1.0)
    coefficient("theta4", 1.0)
    coefficient("u", [(x, y, z, t) -> UX, (x, y, z, t) -> UY])   # velocity vector
    coefficient("gamma", GAM)
    coefficient("epsilon", eps)
    coefficient("delta", 1e-12)     # keeps nhat finite where grad(phis) = 0
    coefficient("ex", [(x, y, z, t) -> 1.0, (x, y, z, t) -> 0.0])   # unit vectors, to
    coefficient("ey", [(x, y, z, t) -> 0.0, (x, y, z, t) -> 1.0])   # pick components

    # The time loop. Every entry is quasi-static; the stepper only sets dt.
    time_interval(T)
    time_stepper(STEP_EULER_IMPLICIT; dt = dt, store_trajectory = snapshots > 0,
                 snapshot_every = max(1, nsteps ÷ max(snapshots, 1)))

    # nhat . grad(v), with nhat from the projected gradient of phis
    ndotgv = "(gx*dot(ex, grad(v)) + gy*dot(ey, grad(v)))*pow(gx*gx + gy*gy + delta, -0.5)"

    # Weak form of one Picard iteration (residual = 0). Terms with phi go to
    # the matrix, the rest to the right-hand side:
    #   theta1*(a0*phi + a1*phio + a2*phioo, v)          time derivative (BDF)
    #   + theta2*(u.grad phi, v)                          advection
    #   + theta3*(gamma*eps*grad phi, grad v)             diffusion
    #   - theta4*(gamma*S(phi, phis)*nhat(phis), grad v)  sharpening, linearized
    # With div(u) = 0 the skew-symmetric form of the derivation is u.grad(phi).
    # Left un-integrated it carries its own outflow: (u.grad phi, v) equals
    # -(u*phi, grad v) + <(u.n)*phi, v> over the whole boundary, the convective
    # outflow flux of the V1 version (inflow rows are Dirichlet, u.n = 0 on top
    # and bottom); the diffusive and sharpening boundary fluxes are dropped
    # ("do-nothing"). V1 took that flux from phis; here it is implicit, since
    # Finch V2 does not allow fields in ds terms.
    weak_form(gx, "gx*v - dot(grad(phis), ex)*v")
    weak_form(gy, "gy*v - dot(grad(phis), ey)*v")
    weak_form(phi, "theta1*(a0*phi + a1*phio + a2*phioo)*v" *
                   " + theta2*dot(u, grad(phi))*v" *
                   " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
                   " - theta4*gamma*" * SHARPENING[form] * "*" * ndotgv)

    first_lag = form == "extrap" ? "2*phio - phioo" : "phio"
    entries = Any[update(phis => first_lag), gx, gy, phi]
    for it = 2:iters
        append!(entries, [update(phis => "phi"), gx, gy, phi])
    end
    append!(entries, [update(phioo => "phio"), update(phio => "phi"),
                      update(a0 => "$(1.5/dt)"), update(a1 => "$(-2.0/dt)"), update(a2 => "$(0.5/dt)"),
                      monitor("mass", "phi")])
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
    y = ctx.grid.nodes[2, :]
    snaps = Vector{Vector{Float64}}()
    times = Float64[]
    mass  = Float64[]
    results = joinpath(proj_dir, "results")
    if ok
        mass = readdlm(joinpath(results, "monitors.csv"), ',', skipstart = 1)[:, 3]
        if snapshots > 0
            _, metas = Finch.read_results_manifest(results)
            snaps = [Finch.load_snapshot_values(results, s)["phi"] for s in metas]
            times = [Float64(s["t"]) for s in metas]
        end
    end
    return (x = x, y = y, phi = phi.values[1, :], dt = dt, nsteps = nsteps,
            snapshots = snaps, snap_times = times, mass = mass, ok = ok)
end
