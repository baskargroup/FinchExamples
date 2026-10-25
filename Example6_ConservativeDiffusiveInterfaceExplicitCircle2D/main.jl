#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, solved explicitly with
# Finch V2 (Julia target).
#
#   theta1*phi_t + theta2*(1/2)*[div(u*phi) + u.grad(phi)]
#     - theta3*div(gamma*eps*grad(phi)) + theta4*div(gamma*(phi - phi^2)*nhat) = 0
#
# phi is the phase field, 0 in one phase and 1 in the other, and nhat is the
# unit interface normal grad(phi)/|grad(phi)|. The diffusion and sharpening
# terms act against each other and hold the interface at a tanh profile of
# fixed width.
#
# Interface normal: recovered at the NODES by an L2 projection of the finite
# element gradient of the previous level, (g_i, v) = (dphi/dx_i, v), one field
# per component (Finch V2 fields are scalar), solved every step before phi;
# nhat = g/sqrt(g.g + delta). Finch V2 has no lumped mass matrix, so the
# projection uses the consistent one (the V1 version of this example lumped it). Taking the normal element by element from
# the gradient of the linear basis gives a field that is constant per element
# and jumps between elements; on this case that normal produced a grid-scale
# zigzag along the leading edge whenever the time integration was accurate
# (RK4, or the semi-implicit schemes) and, on Zalesak's disk (Example 7),
# wiggly slot walls and shed fragments. The nodal normal removes all of it.
#
# Forward Euler in time: every spatial term is evaluated at the previous level
# (Finch V2 does this itself under STEP_EULER_EXPLICIT), so each step is one
# solve with the constant mass matrix.
#
# Test problem: a circle of radius R carried by a uniform velocity. The profile
# translates with almost no change of shape, so
#   phi(x,y,t) = 0.5*(1 - tanh((r - R)/(2*eps)))
# measured from the moving centre is close to the solution. Unlike the 1D case
# it is not exact, because a curved interface is not an equilibrium profile.
#
# Boundaries: the velocity runs straight through the box, in at x = 0 and out
# at x = 1. At the inflow phi = 0 is imposed, since that is what the entering
# fluid carries. At the outflow the weak form keeps the convective boundary
# integral (u.n)*phi*v and drops the diffusive and sharpening ones, the usual
# "do-nothing" condition, so mass leaves at exactly the rate the flow carries
# it out: d/dt(mass) = -oint (u.n) phi ds. That rate is tiny here because
# phi ~ 1e-10 at the wall, but it is real. Dropping the boundary integral
# altogether (NO_BC alone) would instead impose zero total flux and dam the
# flow at the wall; Dirichlet phi = 0 at the OUTFLOW would zero those nodes
# every step and quietly delete that mass. Top and bottom have u.n = 0, so the
# boundary term vanishes there on its own.
#
# Continuous Galerkin, linear elements on a quad mesh.
#==============================================================================#

using Finch
using DelimitedFiles

# results are written to output/, a folder kept in the repo; Finch's generated
# project (and its raw results) goes to output/finch/
const OUT_DIR  = joinpath(@__DIR__, "output")
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

# Physical parameters
const NEL  = 300
const DX   = 1.0/NEL
const UX   = 1.0                # advection velocity components
const UY   = 0.0
const UMAX = sqrt(UX^2 + UY^2)
const GAM  = UMAX               # sharpening parameter gamma (= max |u|)
const EPS  = DX                 # interface thickness (eps/dx >= (u/gamma + 1)/2)
const R    = 0.15               # circle radius
const XC0  = 0.3                # initial circle center
const YC0  = 0.5

# Time stepping. Explicit diffusive limit in 2D (lumped mass):
# dt <= dx^2/(4*gamma*eps). Consistent mass matrix needs roughly 1/3 of that.
const T      = 0.4
const DT     = 0.2*DX*DX/(4*GAM*EPS)
const NSTEPS = Int(round(T/DT))

# The system matrices are mass matrices (symmetric positive definite), so
# plain conjugate gradients converge in a few iterations.
ctx = init_finch("cdi2d"; solver = :cg, solver_tol = 1e-12,
                 target = JuliaCPUTarget(out_dir = PROJ_DIR, func_name = "cdi2d"))

domain(2)
mesh(MESH_QUAD; elsperdim = NEL)   # unit square
function_space(order = 1)

# Boundary ids: 1 = x==0 (left, inflow), 2 = x==1 (right, outflow),
# 3 = y==0 and 4 = y==1 (bottom and top, u.n = 0)

phi = variable("phi")       # the phase field
gx  = variable("gx")        # projected gradient of phi, x component
gy  = variable("gy")        # and y component
test_symbol("v")

# Inflow: the fluid entering at x = 0 carries phi = 0, so that value is imposed.
# The other three walls get no strong condition. What crosses them is set by
# the boundary integral in the weak form below, which acts only at the outflow.
boundary(phi, 1, BC_DIRICHLET, 0.0)     # left, inflow

# Initial condition: equilibrium tanh profile around the circle
initial(phi, (x, y, z, t) -> 0.5*(1 - tanh((sqrt((x - XC0)^2 + (y - YC0)^2) - R)/(2*EPS))))

# Coefficients (theta as in the derivation, dimensional form)
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", [(x, y, z, t) -> UX, (x, y, z, t) -> UY])   # velocity vector
coefficient("gamma", GAM)
coefficient("epsilon", EPS)
coefficient("delta", 1e-12)     # regularizes |grad(phi)| in the bulk
coefficient("ex", [(x, y, z, t) -> 1.0, (x, y, z, t) -> 0.0])   # unit vectors, to
coefficient("ey", [(x, y, z, t) -> 0.0, (x, y, z, t) -> 1.0])   # pick components

# Snapshots for the time series: phi every 40th step, 60 snapshots plus the
# one at t = 0, the 61 the MATLAB scripts read
time_interval(T)
time_stepper(STEP_EULER_EXPLICIT; dt = DT, store_trajectory = true, snapshot_every = 40)

# nhat . grad(v), with nhat from the projected gradient
const NDOTGV = "(gx*dot(ex, grad(v)) + gy*dot(ey, grad(v)))*pow(gx*gx + gy*gy + delta, -0.5)"

# Weak forms (residual = 0). The projection, one per component:
#   (g_i, v) - (dphi/dx_i, v)
# and the phase field, all spatial terms at the previous level:
#   theta1*(Dt phi, v)
#   + theta2*(u.grad phi, v)                          advection
#   + theta3*(gamma*eps*grad phi, grad v)             diffusion
#   - theta4*(gamma*(phi - phi^2)*nhat, grad v)       sharpening, nodal nhat
# With div(u) = 0 the skew-symmetric form of the derivation is u.grad(phi).
# Left un-integrated it carries its own outflow: (u.grad phi, v) equals
# -(u*phi, grad v) + <(u.n)*phi, v> over the whole boundary, which is the
# conservative form plus the convective boundary integral of the V1 version
# (inflow rows are Dirichlet, u.n = 0 on top and bottom). The diffusive and
# sharpening boundary terms are dropped ("do-nothing"). Finch V2 cannot yet
# put the unknown in a ds term of a nonlinear form, so this is also the form
# it can take.
weak_form(gx, "gx*v - dot(grad(phi), ex)*v")
weak_form(gy, "gy*v - dot(grad(phi), ey)*v")
weak_form(phi, "theta1*Dt(phi)*v" *
               " + theta2*dot(u, grad(phi))*v" *
               " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
               " - theta4*gamma*(phi - phi*phi)*" * NDOTGV)

step_schedule(gx, gy, phi, monitor("mass", "phi"))
solve()
Finch.run_emitted_julia!(ctx.target, phi)

# Error check against the translated exact profile
x = ctx.grid.nodes[1, :]
y = ctx.grid.nodes[2, :]
phivals = phi.values[1, :]
xc = XC0 + UX*T
yc = YC0 + UY*T
r  = sqrt.((x .- xc).^2 .+ (y .- yc).^2)
exact = 0.5 .* (1 .- tanh.((r .- R) ./ (2*EPS)))

results = joinpath(PROJ_DIR, "results")
mon = readdlm(joinpath(results, "monitors.csv"), ',', skipstart = 1)
mass0, massT = mon[1, 3], mon[end, 3]

maxerr = maximum(abs.(phivals .- exact))
println("max error = " * string(maxerr))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))
println("relative mass change = " * string((massT - mass0)/mass0))

data = [x y phivals exact]
writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')

# Time series: columns 1-2 = x, y; columns 3.. = phi at each snapshot time
_, snaps = Finch.read_results_manifest(results)
snapshots  = [Finch.load_snapshot_values(results, s)["phi"] for s in snaps]
snap_times = [Float64(s["t"]) for s in snaps]
writedlm(joinpath(OUT_DIR, "solution_timeseries.csv"), hcat(x, y, snapshots...), ',')
writedlm(joinpath(OUT_DIR, "snapshot_times.csv"), snap_times, ',')
println(string(length(snap_times)) * " snapshots written")
