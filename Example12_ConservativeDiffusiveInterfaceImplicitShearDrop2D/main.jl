#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, semi-implicit, with
# Finch V2 (Julia target):
# a drop in a time-reversed shear (vortex) flow.
#
#   theta1*phi_t + theta2*(1/2)*[div(u*phi) + u.grad(phi)]
#     - theta3*div(gamma*eps*grad(phi)) + theta4*div(gamma*(phi - phi^2)*nhat) = 0
#
# The scheme of Examples 10 and 11 on the drop in shear flow test of
# Mirjalili et al., JCP 401 (2020) 109006, Section 4.3.1.
#
# Time discretization (semi-implicit, two Picard iterations per step):
#   - BDF2 for the time derivative, BDF1 on the first step:
#       phi_t ~ (a0*phi + a1*phio + a2*phioo)/dt
#   - advection and diffusion implicit (linear in phi); the velocity is
#     evaluated at the new time level t^{n+1}
#   - sharpening linearized with a lagged field phis, as derived:
#       gamma*(phi - phi^2)*nhat ~ gamma*phi*(1 - phis)*nhat(phis)
#   - iteration 1 uses phis = phio, iteration 2 uses phis = the result of
#     iteration 1. Each iteration is one linear solve (BiCGSTAB); the matrix
#     depends on phis, so it is rebuilt for iteration 2.
#   - one step is a step_schedule, as in Examples 9-11: phi is solved as a
#     quasi-static equation whose time derivative is the BDF formula over the
#     known fields phio and phioo, and update entries move the history, the
#     lagged field and the BDF coefficients along. Under STEP_EULER_IMPLICIT
#     Finch evaluates the solves at t^{n+1}, which sets the velocity's time.
#
# Interface normal: recovered at the NODES by an L2 projection of the finite
# element gradient of phis, one field per component, refreshed before each
# solve (as in Examples 5-11; consistent mass, Finch V2 has no lumped one).
#
# Test problem (Mirjalili et al. 2020, Section 4.3.1): a drop of radius
# R = 0.15 centred at (0.5, 0.75) in the unit square, in the velocity of the
# stream function
#     Psi(x, y, t) = (1/pi) sin^2(pi x) sin^2(pi y) cos(pi t/T),   T = 4,
#     u = -dPsi/dy = -sin^2(pi x) sin(2 pi y) cos(pi t/T)
#     v =  dPsi/dx =  sin(2 pi x) sin^2(pi y) cos(pi t/T)
# The vortex stretches the drop into a thin spiral up to t = T/2, then the
# flow reverses and undoes the stretching. At t = T the exact solution is the
# initial field, so the error is measured against it (the shape error of the
# paper, Eq. 41, is the L1 error below).
#
# Boundaries: no condition on phi and no boundary integral, the zero-total-
# flux (closed box) condition, as in Examples 7 and 11. The normal velocity
# vanishes on the walls, and the discrete mass is conserved.
#
# Velocity: div(u) = 0 (it comes from a stream function), so the div(u)*v
# part of the skew-symmetric advection form is dropped.
#
# Continuous Galerkin, linear elements on a quad mesh.
#==============================================================================#

using Finch
using DelimitedFiles

# Results are written to output/ (Finch's generated project to output/finch/).
# The mesh, the interface parameters and the output folder can be overridden
# from the environment (see paper512.sh); the defaults are the 200 x 200 case.
const OUT_DIR  = joinpath(@__DIR__, get(ENV, "CDI_OUTDIR", "output"))
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

# Physical parameters
const NEL  = parse(Int, get(ENV, "CDI_NEL", "200"))
const DX   = 1.0/NEL
const UMAX = 1.0                # largest |u| over space and time (at t = 0 and T)
const GAM  = parse(Float64, get(ENV, "CDI_GAMMA", "1.0"))*UMAX    # sharpening parameter gamma
const EPS  = parse(Float64, get(ENV, "CDI_EPS_DX", "1.0"))*DX     # interface thickness (eps/dx >= (u/gamma + 1)/2)

# Drop (Mirjalili et al. 2020, Section 4.3.1)
const XC = 0.5
const YC = 0.75
const R  = 0.15
drop_ic(x, y) = 0.5*(1 - tanh((sqrt((x - XC)^2 + (y - YC)^2) - R)/(2*EPS)))

# Time step from the advective CFL number with the largest speed; CFL 0.25
# as in Examples 9-11. One full period of the flow.
const T      = 4.0
const CFL    = 0.25
const NSTEPS = Int(round(T/(CFL*DX/UMAX)))
const DT     = T/NSTEPS             # land exactly on t = T

# The matrix changes every Picard iteration and is not symmetric (advection):
# BiCGSTAB (see Example 10).
ctx = init_finch("cdi2dimplicitsheardrop"; verbose = true, solver = :bicgstab, solver_tol = 1e-12,
                 target = JuliaCPUTarget(out_dir = PROJ_DIR, func_name = "cdi2dimplicitsheardrop"))

domain(2)
mesh(MESH_QUAD; elsperdim = NEL)   # unit square
function_space(order = 1)

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

# Zero total flux through the walls: no boundary condition and no boundary
# integral, so nothing is declared.

# Initial condition, in every level
for f in (phi, phio, phioo, phis)
    initial(f, (x, y, z, t) -> drop_ic(x, y))
end
initial(a0, 1.0/DT); initial(a1, -1.0/DT); initial(a2, 0.0)

# Coefficients (theta as in the derivation, dimensional form)
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
# The vortex, reversed by cos(pi*t/T), at the new time level
coefficient("u", [(x, y, z, t) -> -sin(pi*x)^2*sin(2*pi*y)*cos(pi*t/T),
                  (x, y, z, t) ->  sin(2*pi*x)*sin(pi*y)^2*cos(pi*t/T)])
coefficient("gamma", GAM)
coefficient("epsilon", EPS)
coefficient("delta", 1e-12)     # keeps nhat finite where grad(phis) = 0
coefficient("ex", [(x, y, z, t) -> 1.0, (x, y, z, t) -> 0.0])   # unit vectors, to
coefficient("ey", [(x, y, z, t) -> 0.0, (x, y, z, t) -> 1.0])   # pick components

# The time loop; every entry is quasi-static, the stepper only sets dt.
# Snapshots: 64 over the period plus t = 0; the one at T/2 is the most
# stretched state.
time_interval(T)
time_stepper(STEP_EULER_IMPLICIT; dt = DT, store_trajectory = true,
             snapshot_every = max(1, NSTEPS ÷ 64))

# nhat . grad(v), with nhat from the projected gradient of phis
const NDOTGV = "(gx*dot(ex, grad(v)) + gy*dot(ey, grad(v)))*pow(gx*gx + gy*gy + delta, -0.5)"

# Weak form of one Picard iteration (residual = 0):
#   theta1*(a0*phi + a1*phio + a2*phioo, v)          time derivative (BDF)
#   - theta2/2*(u*phi, grad v)                        advection, conservative half
#   - theta2/2*(phi, u.grad v)                        other half (div(u) = 0)
#   + theta3*(gamma*eps*grad phi, grad v)             diffusion
#   - theta4*(gamma*phi*(1 - phis)*nhat(phis), grad v)  sharpening, lagged
# No boundary term and no boundary condition: zero total flux through the walls.
weak_form(gx, "gx*v - dot(grad(phis), ex)*v")
weak_form(gy, "gy*v - dot(grad(phis), ey)*v")
weak_form(phi, "theta1*(a0*phi + a1*phio + a2*phioo)*v" *
               " - 0.5*theta2*phi*dot(u, grad(v))" *
               " - 0.5*theta2*phi*dot(u, grad(v))" *
               " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
               " - theta4*gamma*phi*(1 - phis)*" * NDOTGV)

# One step: two Picard iterations, then shift the history
step_schedule(update(phis => "phio"), gx, gy, phi,     # iteration 1: lag with the previous level
              update(phis => "phi"),  gx, gy, phi,     # iteration 2: lag with iteration 1
              update(phioo => "phio"), update(phio => "phi"),
              update(a0 => "$(1.5/DT)"), update(a1 => "$(-2.0/DT)"), update(a2 => "$(0.5/DT)"),
              monitor("mass", "phi"))
solve()
Finch.run_emitted_julia!(ctx.target, phi)

# Error after one period, against the initial field
x = ctx.grid.nodes[1, :]
y = ctx.grid.nodes[2, :]
results = joinpath(PROJ_DIR, "results")
_, snaps = Finch.read_results_manifest(results)
snapshots  = [Finch.load_snapshot_values(results, s)["phi"] for s in snaps]
snap_times = [Float64(s["t"]) for s in snaps]
phi0    = snapshots[1]          # the exact solution at t = T
phihalf = snapshots[argmin(abs.(snap_times .- T/2))]   # the most stretched state
phivals = phi.values[1, :]

# Lumped finite element weights (half at an edge node, quarter at a corner)
wts = fill(DX*DX, length(x))
for i in eachindex(x)
    if x[i] < 1e-12 || x[i] > 1 - 1e-12; wts[i] *= 0.5; end
    if y[i] < 1e-12 || y[i] > 1 - 1e-12; wts[i] *= 0.5; end
end
err = phivals .- phi0
maxerr = maximum(abs.(err))
l1err  = sum(wts .* abs.(err))            # shape error, Eq. 41 of the paper
mon = readdlm(joinpath(results, "monitors.csv"), ',', skipstart = 1)
mass0, massT = mon[1, 3], mon[end, 3]
println("mesh " * string(NEL) * " x " * string(NEL) * ", eps/dx = " * string(EPS/DX) * ", gamma/umax = " * string(GAM/UMAX))
println("dt = " * string(DT) * ", steps = " * string(NSTEPS) * " (CFL " * string(CFL) * ", one period T = " * string(T) * ")")
println("max error = " * string(maxerr))
println("L1 (shape) error = " * string(l1err))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))
println("relative mass change = " * string((massT - mass0)/mass0))

data = [x y phivals phi0 phihalf]
writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')

# Time series: columns 1-2 = x, y; columns 3.. = phi at each snapshot time
writedlm(joinpath(OUT_DIR, "solution_timeseries.csv"), hcat(x, y, snapshots...), ',')
writedlm(joinpath(OUT_DIR, "snapshot_times.csv"), snap_times, ',')
println(string(length(snap_times)) * " snapshots written")
