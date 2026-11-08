#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, solved explicitly with
# Finch V2 (branch dev). Run with ./run.sh; the target (CPU Julia, or MLIR on
# the CPU/GPU) is chosen by FINCH_TARGET, see ../common.jl.
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
# element gradient of the previous level, (g, v) = (grad phi, v), a vector
# field solved every step before phi; nhat = g/sqrt(g.g + delta). Finch V2 has
# no lumped mass matrix, so the projection uses the consistent one (the V1
# version lumped it). Taking the normal element by element from the gradient
# of the linear basis gives a field that is constant per element and jumps
# between elements; on the circle that normal produced a grid-scale zigzag
# along the leading edge whenever the time integration was accurate (RK4, or
# the semi-implicit schemes) and, on Zalesak's disk, wiggly slot walls and
# shed fragments. The nodal normal removes all of it (CDI_NORMAL=elementwise
# reproduces the element-wise one).
#
# Forward Euler in time: every spatial term is evaluated at the previous level
# (Finch V2 does this itself under STEP_EULER_EXPLICIT), so each step is one
# solve with the constant mass matrix.
#
# Benchmarks (CDI_CASE):
#   circle     (default) a circle of radius R carried by the uniform velocity
#              u = (1, 0) through the unit square, T = 0.4, mesh 300 x 300.
#              The profile translates with almost no change of shape, so
#              phi = 0.5*(1 - tanh((r - R)/(2*eps))) from the moving centre is
#              close to the solution (not exact: a curved interface is not an
#              equilibrium profile). Boundaries: in at x = 0 (phi = 0 imposed,
#              what the entering fluid carries), out at x = 1, where the weak
#              form keeps the convective boundary integral (u.n)*phi*v and
#              drops the diffusive and sharpening ones (the "do-nothing"
#              condition), so mass leaves exactly at the rate the flow carries
#              it out. Top and bottom have u.n = 0.
#   zalesak    Zalesak's slotted disk in solid-body rotation about the centre
#              of the square, one revolution (T = 2*pi), mesh 200 x 200; the
#              exact solution is the initial field. Closed box: no boundary
#              term and no boundary condition (zero total flux).
#   sheardrop  a drop in the time-reversed vortex of Mirjalili et al. 2020,
#              Section 4.3.1, one period (T = 4), mesh 200 x 200; the exact
#              solution is again the initial field, and the most stretched
#              state (T/2) is saved too. Closed box as for zalesak.
# The closed boxes use the split conservative advection form (div u = 0); the
# circle leaves u.grad(phi) un-integrated so that it carries its own outflow
# flux (Finch V2 cannot put the unknown in a ds term of a nonlinear form).
#
# Environment: CDI_NEL (mesh), CDI_EPS_DX (eps = CDI_EPS_DX*dx, default 1),
# CDI_GAMMA (gamma = CDI_GAMMA*umax; 1.0 when unset), CDI_NORMAL
# (nodal/elementwise), CDI_OUTDIR (default output/<case>/). paper512_*.sh run
# the paper's 512 x 512 parameters, elementwise_normal.sh the comparison.
#
# Continuous Galerkin, linear elements on a quad mesh.
#==============================================================================#

include(joinpath(@__DIR__, "..", "common.jl"))   # Finch, FINCH_TARGET switch, helpers

const CASE = get(ENV, "CDI_CASE", "circle")
CASE in ("circle", "zalesak", "sheardrop") ||
    error("CDI_CASE must be circle, zalesak or sheardrop, got \"$CASE\"")

# Results are written to output/<case>/ (Finch's generated project to its
# finch/ subfolder); CDI_OUTDIR overrides.
const OUT_DIR  = joinpath(@__DIR__, get(ENV, "CDI_OUTDIR", joinpath("output", CASE)))
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

# ---- the benchmark ----------------------------------------------------------------
const NEL  = parse(Int, get(ENV, "CDI_NEL", CASE == "circle" ? "300" : "200"))
const DX   = 1.0/NEL
const T    = CASE == "circle" ? 0.4 : CASE == "zalesak" ? 2*pi : 4.0
const UMAX = CASE == "circle" ? 1.0 : CASE == "zalesak" ? 0.5*sqrt(2) : 1.0   # largest |u|
# sharpening parameter gamma: 1 by default, CDI_GAMMA * umax when given
const GAM  = haskey(ENV, "CDI_GAMMA") ? parse(Float64, ENV["CDI_GAMMA"])*UMAX : 1.0
const EPS  = parse(Float64, get(ENV, "CDI_EPS_DX", "1.0"))*DX   # interface thickness (eps/dx >= (u/gamma + 1)/2)
const NODAL_NORMAL = get(ENV, "CDI_NORMAL", "nodal") != "elementwise"

# circle: centre and radius; it is translated by u*T
const XC0, YC0, R = 0.3, 0.5, 0.15
circle_ic(x, y) = 0.5*(1 - tanh((sqrt((x - XC0)^2 + (y - YC0)^2) - R)/(2*EPS)))

# Zalesak's slotted disk: centre, radius, slot width and the y-range of the
# slot. The slot opens at the bottom of the disk (y = 0.6) and runs up past its
# centre. The slot RECTANGLE used in the distance function is extended 0.05
# below the disk: it removes nothing more from the disk, but if its bottom edge
# sat at the disk's bottom the distance function would see a boundary there and
# the tanh profile would close the slot mouth with a ridge of phi = 0.5.
const ZXC, ZYC, ZR = 0.5, 0.75, 0.15
const SW  = 0.05                # slot width
const SY1 = ZYC - ZR - 0.05     # slot rectangle bottom, below the disk (see above)
const SY2 = ZYC + 0.10          # slot top
# Signed distance to the slotted disk, negative inside. The disk minus the slot
# is "inside disk AND outside slot", whose distance is the larger of the two
# signed distances (exact away from the corners, a good approximation at them).
# Traced by Finch: no data-dependent branches.
function zalesak_sdf(x, y)
    d_disk = sqrt((x - ZXC)^2 + (y - ZYC)^2) - ZR           # negative inside the disk
    ddx = max(abs(x - ZXC) - SW/2, 0.0)
    ddy = max(max(SY1 - y, y - SY2), 0.0)
    outside = sqrt(ddx^2 + ddy^2)
    inside  = min(max(abs(x - ZXC) - SW/2, max(SY1 - y, y - SY2)), 0.0)
    d_slot  = outside + inside                               # negative inside the slot
    return max(d_disk, -d_slot)
end
zalesak_ic(x, y) = 0.5*(1 - tanh(zalesak_sdf(x, y)/(2*EPS)))

# shear drop: the same disk without the slot
const DXC, DYC, DR = 0.5, 0.75, 0.15
drop_ic(x, y) = 0.5*(1 - tanh((sqrt((x - DXC)^2 + (y - DYC)^2) - DR)/(2*EPS)))

const IC = CASE == "circle" ? circle_ic : CASE == "zalesak" ? zalesak_ic : drop_ic
# the velocity: uniform; solid-body rotation; the vortex reversed by cos(pi*t/T)
# (explicit steps read it at the previous level)
const U = CASE == "circle"  ? [(x, y, z, t) -> 1.0, (x, y, z, t) -> 0.0] :
          CASE == "zalesak" ? [(x, y, z, t) -> 0.5 - y, (x, y, z, t) -> x - 0.5] :
                              [(x, y, z, t) -> -sin(pi*x)^2*sin(2*pi*y)*cos(pi*t/T),
                               (x, y, z, t) ->  sin(2*pi*x)*sin(pi*y)^2*cos(pi*t/T)]

# Time stepping. Explicit diffusive limit in 2D (lumped mass):
# dt <= dx^2/(4*gamma*eps). Consistent mass matrix needs roughly 1/3 of that.
# The step count lands exactly on T.
const NSTEPS = Int(round(T/(0.2*DX*DX/(4*GAM*EPS))))
const DT     = T/NSTEPS             # the step Finch takes (Nsteps is what it is given)
# Snapshots for the time series: 60 plus t = 0 (circle: every 40th step; the
# MATLAB scripts read 61), about 60, or 64 over the period with the one at T/2
const SNAP_EVERY = CASE == "circle" ? 40 : CASE == "zalesak" ? max(1, NSTEPS ÷ 60) : max(1, NSTEPS ÷ 64)

# The system matrices are mass matrices (symmetric positive definite): CG. The
# tight tolerance keeps the mass drift of tens of thousands of steps at round-off.
const NAME = CASE == "circle" ? "cdi2d" : CASE == "zalesak" ? "cdi2dzalesak" : "cdi2dsheardrop"
ctx = init_finch(NAME; finch_solver_kwargs(:cg; tol = 1e-12)...,
                 target = finch_target(NAME, PROJ_DIR))

domain(2)
mesh(MESH_QUAD; elsperdim = NEL)   # unit square
function_space(order = 1)

# Boundary ids: 1 = x==0 (left), 2 = x==1 (right), 3 = y==0 and 4 = y==1

phi = variable("phi")       # the phase field
g   = variable("g"; kind = VAR_VECTOR, components = 2)   # projected gradient of phi (nodal normal)
test_symbol("v")

# circle: inflow at x = 0 carries phi = 0; the other walls get no strong
# condition (the outflow flux is in the weak form). Closed boxes: nothing.
CASE == "circle" && boundary(phi, 1, BC_DIRICHLET, 0.0)

# Initial condition: equilibrium tanh profile across the shape's boundary
initial(phi, (x, y, z, t) -> IC(x, y))

# Coefficients (theta as in the derivation, dimensional form)
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", U)             # velocity vector
coefficient("gamma", GAM)
coefficient("epsilon", EPS)
coefficient("delta", 1e-12)     # regularizes |grad(phi)| in the bulk

time_interval(T)
time_stepper(STEP_EULER_EXPLICIT; Nsteps = NSTEPS, store_trajectory = true, snapshot_every = SNAP_EVERY)

# nhat . grad(v): from the projected gradient (nodal), or element by element
# from the gradient of the linear basis (CDI_NORMAL=elementwise); common.jl
const NDOTGV = NODAL_NORMAL ? NDOTGV_NODAL_2D : NDOTGV_ELEMENT_2D

# Weak forms (residual = 0). The projection (vector unknown g, vector test v):
#   (g, v) - (grad phi, v)
# and the phase field, all spatial terms at the previous level:
#   theta1*(Dt phi, v)
#   advection: circle      + theta2*(u.grad phi, v)     un-integrated, carries the outflow:
#                            (u.grad phi, v) = -(u*phi, grad v) + <(u.n)*phi, v>
#              closed box  - theta2/2*(u*phi, grad v) - theta2/2*(phi, u.grad v)
#   + theta3*(gamma*eps*grad phi, grad v)             diffusion
#   - theta4*(gamma*(phi - phi^2)*nhat, grad v)       sharpening
const ADVECTION = CASE == "circle" ? " + theta2*dot(u, grad(phi))*v" :
                  " - 0.5*theta2*phi*dot(u, grad(v)) - 0.5*theta2*phi*dot(u, grad(v))"
weak_form(g, "dot(g, v) - dot(grad(phi), v)")
weak_form(phi, "theta1*Dt(phi)*v" * ADVECTION *
               " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
               " - theta4*gamma*(phi - phi*phi)*" * NDOTGV)

if CASE == "zalesak"
    println("normal: " * (NODAL_NORMAL ? "nodal" : "element-wise"))
    println("dt = " * string(DT) * ", steps = " * string(NSTEPS) * " (one revolution, T = 2*pi)")
elseif CASE == "sheardrop"
    println("mesh " * string(NEL) * " x " * string(NEL) * ", eps/dx = " * string(EPS/DX) * ", gamma/umax = " * string(GAM/UMAX))
    println("dt = " * string(DT) * ", steps = " * string(NSTEPS) * " (one period, T = " * string(T) * ")")
else
    println("dt = $DT, steps = $NSTEPS")
end
if NODAL_NORMAL
    step_schedule(g, phi, monitor("mass", "phi"))
else
    step_schedule(phi, monitor("mass", "phi"))
end
solve()
finch_run!(ctx, phi)

# ---- errors and output: the final state is the last snapshot ----------------------
x = ctx.grid.nodes[1, :]
y = ctx.grid.nodes[2, :]
snap_times, fields = read_snapshots(PROJ_DIR)
check_final_time(snap_times, T)
phivals = final_values(fields, "phi")
mass0, massT, _ = mass_change(joinpath(PROJ_DIR, "results"))

if CASE == "circle"
    # against the translated profile
    xc, yc = XC0 + 1.0*T, YC0 + 0.0*T
    r  = sqrt.((x .- xc).^2 .+ (y .- yc).^2)
    exact = 0.5 .* (1 .- tanh.((r .- R) ./ (2*EPS)))
    println("max error = " * string(maximum(abs.(phivals .- exact))))
    data = [x y phivals exact]
else
    # after one revolution or one period, against the initial field; L1 with
    # lumped weights (half at an edge node, quarter at a corner)
    phi0 = fields["phi"][:, 1]
    wts  = lumped_weights(x, y, DX)
    err  = phivals .- phi0
    println("max error = " * string(maximum(abs.(err))))
    if CASE == "zalesak"
        println("L1 error = " * string(sum(wts .* abs.(err))))          # the usual Zalesak measure
        data = [x y phivals phi0]
    else
        println("L1 (shape) error = " * string(sum(wts .* abs.(err))))  # shape error, Eq. 41 of the paper
        phihalf = fields["phi"][:, argmin(abs.(snap_times .- T/2))]     # the most stretched state
        data = [x y phivals phi0 phihalf]
    end
end
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))
println("relative mass change = " * string((massT - mass0)/mass0))

writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')
# Time series: columns 1-2 = x, y; columns 3.. = phi at each snapshot time
write_timeseries(OUT_DIR, "solution", [x y], fields["phi"])
write_snapshot_times(OUT_DIR, snap_times)
println(string(length(snap_times)) * " snapshots written")
