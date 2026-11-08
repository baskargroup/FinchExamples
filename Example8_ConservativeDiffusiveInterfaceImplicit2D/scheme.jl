#==============================================================================#
# The semi-implicit 2D CDI scheme as a function, so that main.jl and
# sharpening_sweep.jl run the same code. Finch V2 builds a fresh context (and
# generated project) for each call.
#
#   run_case(CFL, iters; case, form, nel, proj_dir, snapshots)
#
# case selects the benchmark (the flow, the initial shape, the walls, T):
#   "circle"     a circle carried by u = (1, 0) through the square, inflow on
#                the left, outflow on the right (the default; mesh 300 x 300)
#   "zalesak"    Zalesak's slotted disk in solid-body rotation, one revolution,
#                closed box (mesh 200 x 200)
#   "sheardrop"  a drop in the time-reversed vortex of Mirjalili et al. 2020,
#                Section 4.3.1, one period, closed box (mesh 200 x 200)
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
# Returns (x, y, phi(T), dt, nsteps, snapshots, snap_times, mass, eps, gam,
# umax, ok); mass is int phi after every step (index 1 is t = 0).
#
# Finch V2 (branch dev): the level two steps back is Finch's own history
# `old(phi, 2)` and the BDF coefficients are coefficients of t (BDF1 on the
# first step, BDF2 after), so only phio and phis remain as fields. The target
# (CPU Julia, or MLIR on the CPU/GPU) is chosen by FINCH_TARGET, see ../common.jl.
#
# Environment: CDI_GAMMA (gamma = CDI_GAMMA * umax; 1.0 when unset) and
# CDI_EPS_DX (eps = CDI_EPS_DX * dx, default 1.0), as in the explicit example.
#==============================================================================#

include(joinpath(@__DIR__, "..", "common.jl"))   # Finch, FINCH_TARGET switch, helpers

const CASES = ("circle", "zalesak", "sheardrop")

const SHARPENING = Dict(
    "derived"  => "phi*(1 - phis)",
    "lagged"   => "phis*(1 - phis)",
    "reversed" => "phis*(1 - phi)",
    "newton"   => "((1 - 2*phis)*phi + phis*phis)",
    "extrap"   => "phis*(1 - phis)")

# ---- the three benchmarks ------------------------------------------------------
# circle: centre and radius; the circle is translated by u*T
const XC0, YC0, R = 0.3, 0.5, 0.15
circle_ic(x, y, eps) = 0.5*(1 - tanh((sqrt((x - XC0)^2 + (y - YC0)^2) - R)/(2*eps)))

# Zalesak's slotted disk: the slot rectangle extends 0.05 below the disk so that
# the slot mouth is open (see the explicit example)
const ZXC, ZYC, ZR = 0.5, 0.75, 0.15
const SW  = 0.05                # slot width
const SY1 = ZYC - ZR - 0.05     # slot rectangle bottom, below the disk
const SY2 = ZYC + 0.10          # slot top
function zalesak_sdf(x, y)      # signed distance, negative inside (traced by Finch)
    d_disk = sqrt((x - ZXC)^2 + (y - ZYC)^2) - ZR
    ddx = max(abs(x - ZXC) - SW/2, 0.0)
    ddy = max(max(SY1 - y, y - SY2), 0.0)
    outside = sqrt(ddx^2 + ddy^2)
    inside  = min(max(abs(x - ZXC) - SW/2, max(SY1 - y, y - SY2)), 0.0)
    d_slot  = outside + inside
    return max(d_disk, -d_slot)
end
zalesak_ic(x, y, eps) = 0.5*(1 - tanh(zalesak_sdf(x, y)/(2*eps)))

# shear drop: the same disk without the slot, in the vortex reversed at T/2
const DXC, DYC, DR = 0.5, 0.75, 0.15
drop_ic(x, y, eps) = 0.5*(1 - tanh((sqrt((x - DXC)^2 + (y - DYC)^2) - DR)/(2*eps)))

# Per case: default mesh, final time, largest speed (for gamma and the CFL
# step), how nsteps is rounded, the snapshot count, the velocity closures, the
# initial condition, whether the left wall is an inflow, and the advection form.
# The circle leaves u.grad(phi) un-integrated so that it carries its own
# outflow flux; the closed boxes use the split conservative form (div u = 0)
# with no boundary term and no boundary condition.
function case_setup(case::String, T::Float64)
    if case == "circle"
        return (nel = 300, T = 0.4, umax = 1.0, round_steps = false, snapshots = 60,
                u = [(x, y, z, t) -> 1.0, (x, y, z, t) -> 0.0],
                ic = circle_ic, inflow = true,
                advection = " + theta2*dot(u, grad(phi))*v")
    elseif case == "zalesak"
        return (nel = 200, T = 2*pi, umax = 0.5*sqrt(2), round_steps = true, snapshots = 60,
                u = [(x, y, z, t) -> 0.5 - y, (x, y, z, t) -> x - 0.5],
                ic = zalesak_ic, inflow = false,
                advection = " - 0.5*theta2*phi*dot(u, grad(v)) - 0.5*theta2*phi*dot(u, grad(v))")
    elseif case == "sheardrop"
        return (nel = 200, T = 4.0, umax = 1.0, round_steps = true, snapshots = 64,
                u = [(x, y, z, t) -> -sin(pi*x)^2*sin(2*pi*y)*cos(pi*t/T),
                     (x, y, z, t) ->  sin(2*pi*x)*sin(pi*y)^2*cos(pi*t/T)],
                ic = drop_ic, inflow = false,
                advection = " - 0.5*theta2*phi*dot(u, grad(v)) - 0.5*theta2*phi*dot(u, grad(v))")
    end
    error("CDI_CASE must be circle, zalesak or sheardrop, got \"$case\"")
end
case_T(case::String) = case == "circle" ? 0.4 : case == "zalesak" ? 2*pi : case == "sheardrop" ? 4.0 :
                       error("CDI_CASE must be circle, zalesak or sheardrop, got \"$case\"")
case_nel(case::String) = case_setup(case, case_T(case)).nel

function run_case(CFL, iters; case = "circle", form = "derived", nel = case_nel(case),
                  proj_dir = mktempdir(), snapshots = -1, verbose = false)
    T  = case_T(case)
    cs = case_setup(case, T)
    snapshots < 0 && (snapshots = cs.snapshots)
    dx   = 1.0/nel
    eps  = parse(Float64, get(ENV, "CDI_EPS_DX", "1.0"))*dx      # interface thickness (eps/dx >= (u/gamma + 1)/2)
    gam  = haskey(ENV, "CDI_GAMMA") ? parse(Float64, ENV["CDI_GAMMA"])*cs.umax : 1.0   # sharpening parameter
    # dt from the advective CFL number with the largest speed; the circle rounds
    # dt down so that T is reached exactly, the closed boxes round to nearest
    nsteps = cs.round_steps ? Int(round(T/(CFL*dx/cs.umax))) : ceil(Int, T/(CFL*dx/cs.umax) - 1e-9)
    dt = T/nsteps
    mkpath(proj_dir)

    # The matrix changes every Picard iteration and is not symmetric
    # (advection): BiCGSTAB. A loose tolerance would show up directly as a
    # mass drift of that size.
    ctx = init_finch("cdi2dimplicit"; verbose = verbose, finch_solver_kwargs(:bicgstab; tol = 1e-12)...,
                     target = finch_target("cdi2dimplicit", proj_dir))
    domain(2)
    mesh(MESH_QUAD; elsperdim = nel)   # unit square
    function_space(order = 1)

    # Boundary ids: 1 = x==0 (left), 2 = x==1 (right), 3 = y==0 and 4 = y==1

    phi   = variable("phi")     # unknown: the current Picard iterate; old(phi, 2) is two levels back
    phio  = variable("phio")    # previous time level (old(phi) would be the Picard iterate)
    phis  = variable("phis")    # lagged field for the sharpening term
    g     = variable("g"; kind = VAR_VECTOR, components = 2)   # projected gradient of phis (nodal normal)
    test_symbol("v")

    # Inflow (circle only): the fluid entering at x = 0 carries phi = 0. The
    # closed boxes have zero total flux through the walls: nothing is declared.
    cs.inflow && boundary(phi, 1, BC_DIRICHLET, 0.0)

    # Initial condition: equilibrium tanh profile across the shape, in every level
    ic = cs.ic
    for f in (phi, phio, phis)
        initial(f, (x, y, z, t) -> ic(x, y, eps))
    end
    # BDF coefficients over dt: BDF1 on the first step (the solves of step n see
    # t = n*dt), BDF2 after
    coefficient("a0", (x, y, z, t) -> ifelse(t < 1.5*dt,  1.0/dt,  1.5/dt))
    coefficient("a1", (x, y, z, t) -> ifelse(t < 1.5*dt, -1.0/dt, -2.0/dt))
    coefficient("a2", (x, y, z, t) -> ifelse(t < 1.5*dt,  0.0,     0.5/dt))

    # Coefficients (theta as in the derivation, dimensional form)
    coefficient("theta1", 1.0)
    coefficient("theta2", 1.0)
    coefficient("theta3", 1.0)
    coefficient("theta4", 1.0)
    coefficient("u", cs.u)          # velocity vector
    coefficient("gamma", gam)
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
    #   advection (per case, see case_setup)
    #   + theta3*(gamma*eps*grad phi, grad v)             diffusion
    #   - theta4*(gamma*S(phi, phis)*nhat(phis), grad v)  sharpening, linearized
    # The projection of the lagged gradient (vector unknown g, vector test v):
    #   (g, v) - (grad phis, v)
    weak_form(g, "dot(g, v) - dot(grad(phis), v)")
    weak_form(phi, "theta1*(a0*phi + a1*phio + a2*old(phi, 2))*v" *
                   cs.advection *
                   " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
                   " - theta4*gamma*" * SHARPENING[form] * "*" * NDOTGV_NODAL_2D)

    first_lag = form == "extrap" ? "2*phio - old(phi, 2)" : "phio"
    entries = Any[update(phis => first_lag), g, phi]
    for it = 2:iters
        append!(entries, [update(phis => "phi"), g, phi])
    end
    append!(entries, [update(phio => "phi"), monitor("mass", "phi")])
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
    y = ctx.grid.nodes[2, :]
    snaps = Vector{Vector{Float64}}()
    times = Float64[]
    mass  = Float64[]
    phiT  = fill(NaN, length(x))
    if ok
        mass = readdlm(joinpath(proj_dir, "results", "monitors.csv"), ',', skipstart = 1)[:, 3]
        # the final state is the last entry of results/ (both targets write it)
        times, fields = read_snapshots(proj_dir)
        check_final_time(times, T)
        phiT = final_values(fields, "phi")
        snaps = snapshots > 0 ? [fields["phi"][:, k] for k in eachindex(times)] : snaps
        snapshots > 0 || (times = Float64[])
    end
    return (x = x, y = y, phi = phiT, dt = dt, nsteps = nsteps,
            snapshots = snaps, snap_times = times, mass = mass,
            eps = eps, gam = gam, umax = cs.umax, T = T, ok = ok)
end
