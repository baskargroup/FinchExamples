#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, solved explicitly with
# Finch V2 (Julia target):
# Zalesak's slotted disk in solid-body rotation.
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
# per component, solved every step before phi (Examples 5 and 6 do the same).
# Finch V2 has no lumped mass matrix, so the projection uses the consistent
# one (the V1 version of this example lumped it). Reference Examples 10-16
# take the normal element by element from the gradient of the linear basis;
# that normal is constant inside each element and jumps between elements, and
# on this case it produced grid-scale wiggles along the slot walls, shed
# fragments and a widened slot (L1 error 1.2e-2, slot 0.059 wide instead of
# 0.050; elementwise_normal.sh reproduces it, CDI_NORMAL=elementwise). With the nodal normal the
# walls are smooth, nothing is shed and the slot keeps its width (L1 error
# 3.4e-3). Mirjalili et al. (JCP 2020) likewise use a continuous normal, from
# central differences at cell centres.
#
# Scheme: the same as Example 6. Forward Euler in time, so each step is one
# solve with the constant mass matrix; every spatial term is evaluated at the
# previous level (Finch V2 does this itself under STEP_EULER_EXPLICIT).
#
# Test problem: a disk of radius R with a rectangular slot cut into it, carried
# by the rotation u = (0.5 - y, x - 0.5) about the centre of the unit square.
# One revolution takes T = 2*pi, after which the exact solution is the initial
# field, so the error is measured against the initial condition. The disk never
# comes within 0.1 of a wall. The corners of the slot are not an equilibrium
# shape for a diffuse interface, so some rounding over a distance ~eps is
# expected; the body of the disk should return to where it started.
#
# Boundaries: the rotation crosses all four walls, but the disk stays 0.1 away
# from them and phi is ~1e-9 there. No condition is imposed on phi and no
# boundary integral is kept, which is the zero-total-flux (closed box) condition:
# the discrete mass is then conserved exactly. (Dirichlet phi = 0 instead would
# delete the tanh tails that touch the wall: measured -5e-3 of the mass on a
# 50 x 50 mesh and -1e-5 on 100 x 100, where the tails are longer.)
#
# Velocity: div(u) = 0 for solid-body rotation, so the div(u)*v part of the
# skew-symmetric form is dropped, and its two halves, each integrated by
# parts, are both -(1/2)*(u*phi, grad v): together the conservative form,
# with no boundary integral (closed box).
#
# Continuous Galerkin, linear elements on a quad mesh.
#==============================================================================#

using Finch
using DelimitedFiles

# Results are written to output/ (Finch's generated project to output/finch/).
# The mesh, the interface parameters, the normal and the output folder can be
# overridden from the environment (see paper512.sh and elementwise_normal.sh);
# the defaults are the 200 x 200 case with the nodal normal.
const OUT_DIR  = joinpath(@__DIR__, get(ENV, "CDI_OUTDIR", "output"))
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

const NEL  = parse(Int, get(ENV, "CDI_NEL", "200"))
const DX   = 1.0/NEL
const UMAX = 0.5*sqrt(2)        # largest |u| in the square (at the corners)
# sharpening parameter gamma: 1 by default, CDI_GAMMA * umax when given
const GAM  = haskey(ENV, "CDI_GAMMA") ? parse(Float64, ENV["CDI_GAMMA"])*UMAX : 1.0
const EPS  = parse(Float64, get(ENV, "CDI_EPS_DX", "1.0"))*DX   # interface thickness (eps/dx >= (u/gamma + 1)/2)
const NODAL_NORMAL = get(ENV, "CDI_NORMAL", "nodal") != "elementwise"

# Slotted disk: centre, radius, slot width and the y-range of the slot. The
# slot opens at the bottom of the disk (y = 0.6) and runs up past its centre.
# The slot RECTANGLE used in the distance function is extended 0.05 below the
# disk: it removes nothing more from the disk, but if its bottom edge sat at
# the disk's bottom the distance function would see a boundary there and the
# tanh profile would close the slot mouth with a ridge of phi = 0.5.
const XC  = 0.5
const YC  = 0.75
const R   = 0.15
const SW  = 0.05                # slot width
const SY1 = YC - R - 0.05       # slot rectangle bottom, below the disk (see above)
const SY2 = YC + 0.10           # slot top

# Signed distance to the slotted disk, negative inside. The disk minus the slot
# is "inside disk AND outside slot", whose distance is the larger of the two
# signed distances (exact away from the corners, a good approximation at them).
function zalesak_sdf(x, y)
    d_disk = sqrt((x - XC)^2 + (y - YC)^2) - R           # negative inside the disk
    # signed distance to the slot rectangle [xc-sw/2, xc+sw/2] x [sy1, sy2]
    ddx = max(abs(x - XC) - SW/2, 0.0)
    ddy = max(max(SY1 - y, y - SY2), 0.0)
    outside = sqrt(ddx^2 + ddy^2)
    inside  = min(max(abs(x - XC) - SW/2, max(SY1 - y, y - SY2)), 0.0)
    d_slot  = outside + inside                            # negative inside the slot
    return max(d_disk, -d_slot)
end

# Initial condition: equilibrium tanh profile across the shape's boundary
zalesak_ic(x, y) = 0.5*(1 - tanh(zalesak_sdf(x, y)/(2*EPS)))

# Time stepping. Explicit diffusive limit in 2D (lumped mass):
# dt <= dx^2/(4*gamma*eps). Consistent mass matrix needs roughly 1/3 of that.
# One full revolution.
const T      = 2*pi
const NSTEPS = Int(round(T/(0.2*DX*DX/(4*GAM*EPS))))
const DT     = T/NSTEPS             # land exactly on one revolution

# The system matrices are mass matrices (symmetric positive definite): CG. The
# tight tolerance keeps the mass drift of 25,000 steps at round-off.
ctx = init_finch("cdi2dzalesak"; solver = :cg, solver_tol = 1e-12,
                 target = JuliaCPUTarget(out_dir = PROJ_DIR, func_name = "cdi2dzalesak"))

domain(2)
mesh(MESH_QUAD; elsperdim = NEL)   # unit square
function_space(order = 1)

phi = variable("phi")       # the phase field
gx  = variable("gx")        # projected gradient of phi, x component
gy  = variable("gy")        # and y component
test_symbol("v")

# Zero total flux through the walls: no boundary condition and no boundary
# integral, so nothing is declared.

# Initial condition, traced by Finch (no data-dependent branches above)
initial(phi, (x, y, z, t) -> zalesak_ic(x, y))

# Coefficients (theta as in the derivation, dimensional form)
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", [(x, y, z, t) -> 0.5 - y, (x, y, z, t) -> x - 0.5])   # solid-body rotation
coefficient("gamma", GAM)
coefficient("epsilon", EPS)
coefficient("delta", 1e-12)     # regularizes |grad(phi)| in the bulk
coefficient("ex", [(x, y, z, t) -> 1.0, (x, y, z, t) -> 0.0])   # unit vectors, to
coefficient("ey", [(x, y, z, t) -> 0.0, (x, y, z, t) -> 1.0])   # pick components

# Snapshots for the time series: about 60, plus the one at t = 0
time_interval(T)
time_stepper(STEP_EULER_EXPLICIT; dt = DT, store_trajectory = true,
             snapshot_every = max(1, NSTEPS ÷ 60))

# nhat . grad(v): from the projected gradient (nodal), or element by element
# from the gradient of the linear basis (CDI_NORMAL=elementwise)
const NDOTGV = NODAL_NORMAL ?
    "(gx*dot(ex, grad(v)) + gy*dot(ey, grad(v)))*pow(gx*gx + gy*gy + delta, -0.5)" :
    "dot(grad(phi), grad(v))*pow(dot(grad(phi), grad(phi)) + delta, -0.5)"

# Weak forms (residual = 0). The projection, one per component:
#   (g_i, v) - (dphi/dx_i, v)
# and the phase field, all spatial terms at the previous level:
#   theta1*(Dt phi, v)
#   - theta2/2*(u*phi, grad v)                        conservative half
#   - theta2/2*(phi, u.grad v)                        other half (div(u) = 0)
#   + theta3*(gamma*eps*grad phi, grad v)             diffusion
#   - theta4*(gamma*(phi - phi^2)*nhat, grad v)       sharpening
# No boundary term and no boundary condition: zero total flux through the walls.
weak_form(gx, "gx*v - dot(grad(phi), ex)*v")
weak_form(gy, "gy*v - dot(grad(phi), ey)*v")
weak_form(phi, "theta1*Dt(phi)*v" *
               " - 0.5*theta2*phi*dot(u, grad(v))" *
               " - 0.5*theta2*phi*dot(u, grad(v))" *
               " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
               " - theta4*gamma*(phi - phi*phi)*" * NDOTGV)

if NODAL_NORMAL
    step_schedule(gx, gy, phi, monitor("mass", "phi"))
else
    step_schedule(phi, monitor("mass", "phi"))
end
solve()
Finch.run_emitted_julia!(ctx.target, phi)

# Error after one revolution, against the initial field
x = ctx.grid.nodes[1, :]
y = ctx.grid.nodes[2, :]
results = joinpath(PROJ_DIR, "results")
_, snaps = Finch.read_results_manifest(results)
snapshots  = [Finch.load_snapshot_values(results, s)["phi"] for s in snaps]
snap_times = [Float64(s["t"]) for s in snaps]
phi0    = snapshots[1]          # the exact solution after one revolution
phivals = phi.values[1, :]

# Lumped finite element weights (half at an edge node, quarter at a corner)
wts = fill(DX*DX, length(x))
for i in eachindex(x)
    if x[i] < 1e-12 || x[i] > 1 - 1e-12; wts[i] *= 0.5; end
    if y[i] < 1e-12 || y[i] > 1 - 1e-12; wts[i] *= 0.5; end
end
err = phivals .- phi0
maxerr = maximum(abs.(err))
l1err  = sum(wts .* abs.(err))            # the usual Zalesak measure
mon = readdlm(joinpath(results, "monitors.csv"), ',', skipstart = 1)
mass0, massT = mon[1, 3], mon[end, 3]

println("normal: " * (NODAL_NORMAL ? "nodal" : "element-wise"))
println("dt = " * string(DT) * ", steps = " * string(NSTEPS) * " (one revolution, T = 2*pi)")
println("max error = " * string(maxerr))
println("L1 error = " * string(l1err))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))
println("relative mass change = " * string((massT - mass0)/mass0))

data = [x y phivals phi0]
writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')

# Time series: columns 1-2 = x, y; columns 3.. = phi at each snapshot time
writedlm(joinpath(OUT_DIR, "solution_timeseries.csv"), hcat(x, y, snapshots...), ',')
writedlm(joinpath(OUT_DIR, "snapshot_times.csv"), snap_times, ',')
println(string(length(snap_times)) * " snapshots written")
