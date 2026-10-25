#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, solved explicitly with
# Finch V2 (Julia target):
# a drop in a time-reversed shear (vortex) flow.
#
#   theta1*phi_t + theta2*(1/2)*[div(u*phi) + u.grad(phi)]
#     - theta3*div(gamma*eps*grad(phi)) + theta4*div(gamma*(phi - phi^2)*nhat) = 0
#
# phi is the phase field, 0 in one phase and 1 in the other, and nhat is the
# unit interface normal grad(phi)/|grad(phi)|. The diffusion and sharpening
# terms act against each other and hold the interface at a tanh profile of
# fixed width.
#
# Scheme: the same as Examples 6 and 7. Forward Euler in time, so each step is
# one solve with the constant mass matrix; every spatial term, and the
# velocity, is evaluated at the previous level (Finch V2 does this itself
# under STEP_EULER_EXPLICIT).
#
# Interface normal: recovered at the NODES by an L2 projection of the
# finite element gradient of the previous level, solved every step before phi
# (as in Examples 5-7; consistent mass, Finch V2 has no lumped one).
#
# Test problem (Mirjalili et al., JCP 401 (2020) 109006, Section 4.3.1): a
# drop of radius R = 0.15 centred at (0.5, 0.75) in the unit square, in the
# velocity of the stream function
#     Psi(x, y, t) = (1/pi) sin^2(pi x) sin^2(pi y) cos(pi t/T),   T = 4,
#     u = -dPsi/dy = -sin^2(pi x) sin(2 pi y) cos(pi t/T)
#     v =  dPsi/dx =  sin(2 pi x) sin^2(pi y) cos(pi t/T)
# The vortex stretches the drop into a thin spiral up to t = T/2, then the
# flow reverses and undoes the stretching. At t = T the exact solution is the
# initial field, so the error is measured against it (the shape error of the
# paper, Eq. 41, is the L1 error below). Example 12 solves the same problem
# with the semi-implicit scheme.
#
# Boundaries: no condition on phi and no boundary integral, the zero-total-
# flux (closed box) condition, as in Example 7. The normal velocity vanishes
# on the walls, and the discrete mass is conserved.
#
# Velocity: div(u) = 0 (it comes from a stream function), so the div(u)*v
# part of the skew-symmetric form is dropped.
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

const NEL  = parse(Int, get(ENV, "CDI_NEL", "200"))
const DX   = 1.0/NEL
const UMAX = 1.0                # largest |u| over space and time (at t = 0 and T)
const GAM  = parse(Float64, get(ENV, "CDI_GAMMA", "1.0"))*UMAX    # sharpening parameter gamma
const EPS  = parse(Float64, get(ENV, "CDI_EPS_DX", "1.0"))*DX     # interface thickness (eps/dx >= (u/gamma + 1)/2)

# Drop (Mirjalili et al. 2020, Section 4.3.1)
const XC = 0.5
const YC = 0.75
const R  = 0.15

# Time stepping. Explicit diffusive limit in 2D (lumped mass):
# dt <= dx^2/(4*gamma*eps). Consistent mass matrix needs roughly 1/3 of that.
# One full period of the flow.
const T      = 4.0
const NSTEPS = Int(round(T/(0.2*DX*DX/(4*GAM*EPS))))
const DT     = T/NSTEPS             # land exactly on t = T

# The system matrices are mass matrices (symmetric positive definite): CG. The
# tight tolerance keeps the mass drift of 16,000 steps at round-off.
ctx = init_finch("cdi2dsheardrop"; solver = :cg, solver_tol = 1e-12,
                 target = JuliaCPUTarget(out_dir = PROJ_DIR, func_name = "cdi2dsheardrop"))

domain(2)
mesh(MESH_QUAD; elsperdim = NEL)   # unit square
function_space(order = 1)

phi = variable("phi")       # the phase field
gx  = variable("gx")        # projected gradient of phi, x component
gy  = variable("gy")        # and y component
test_symbol("v")

# Zero total flux through the walls: no boundary condition and no boundary
# integral, so nothing is declared.

# Initial condition: equilibrium tanh profile across the drop boundary
initial(phi, (x, y, z, t) -> 0.5*(1 - tanh((sqrt((x - XC)^2 + (y - YC)^2) - R)/(2*EPS))))

# Coefficients (theta as in the derivation, dimensional form)
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
# The vortex, reversed by cos(pi*t/T); explicit steps read it at the previous level
coefficient("u", [(x, y, z, t) -> -sin(pi*x)^2*sin(2*pi*y)*cos(pi*t/T),
                  (x, y, z, t) ->  sin(2*pi*x)*sin(pi*y)^2*cos(pi*t/T)])
coefficient("gamma", GAM)
coefficient("epsilon", EPS)
coefficient("delta", 1e-12)     # regularizes |grad(phi)| in the bulk
coefficient("ex", [(x, y, z, t) -> 1.0, (x, y, z, t) -> 0.0])   # unit vectors, to
coefficient("ey", [(x, y, z, t) -> 0.0, (x, y, z, t) -> 1.0])   # pick components

# Snapshots for the time series: 64 over the period plus t = 0; the one at
# T/2 is the most stretched state
const SNAP_EVERY = max(1, NSTEPS ÷ 64)
time_interval(T)
time_stepper(STEP_EULER_EXPLICIT; dt = DT, store_trajectory = true, snapshot_every = SNAP_EVERY)

# nhat . grad(v), with nhat from the projected gradient
const NDOTGV = "(gx*dot(ex, grad(v)) + gy*dot(ey, grad(v)))*pow(gx*gx + gy*gy + delta, -0.5)"

# Weak forms (residual = 0). The projection, one per component:
#   (g_i, v) - (dphi/dx_i, v)
# and the phase field, all spatial terms at the previous level:
#   theta1*(Dt phi, v)
#   - theta2/2*(u*phi, grad v)                        conservative half
#   - theta2/2*(phi, u.grad v)                        other half (div(u) = 0)
#   + theta3*(gamma*eps*grad phi, grad v)             diffusion
#   - theta4*(gamma*(phi - phi^2)*nhat, grad v)       sharpening, nodal nhat
weak_form(gx, "gx*v - dot(grad(phi), ex)*v")
weak_form(gy, "gy*v - dot(grad(phi), ey)*v")
weak_form(phi, "theta1*Dt(phi)*v" *
               " - 0.5*theta2*phi*dot(u, grad(v))" *
               " - 0.5*theta2*phi*dot(u, grad(v))" *
               " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
               " - theta4*gamma*(phi - phi*phi)*" * NDOTGV)

step_schedule(gx, gy, phi, monitor("mass", "phi"))
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
println("dt = " * string(DT) * ", steps = " * string(NSTEPS) * " (one period, T = " * string(T) * ")")
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
