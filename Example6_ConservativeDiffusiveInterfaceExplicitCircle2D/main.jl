#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, solved explicitly in Finch.
#
#   theta1*phi_t + theta2*(1/2)*[div(u*phi) + u.grad(phi)]
#     - theta3*div(gamma*eps*grad(phi)) + theta4*div(gamma*(phi - phi^2)*nhat) = 0
#
# phi is the phase field, 0 in one phase and 1 in the other, and nhat is the
# unit interface normal grad(phi)/|grad(phi)|. The diffusion and sharpening
# terms act against each other and hold the interface at a tanh profile of
# fixed width.
#
# Interface normal: recovered at the NODES by a lumped L2 projection of the
# finite element gradient of the previous level, stored in a vector variable
# and refreshed after every step. Taking the normal element by element from
# the gradient of the linear basis gives a field that is constant per element
# and jumps between elements; on this case that normal produced a grid-scale
# zigzag along the leading edge whenever the time integration was accurate
# (RK4, or the semi-implicit schemes) and, on Zalesak's disk (Example 7),
# wiggly slot walls and shed fragments. The nodal normal removes all of it.
#
# Forward Euler in time. The unknown phi then appears only in the time term,
# so each step is one linear solve with the constant mass matrix.
# Every spatial term is evaluated at the previous level, held in phio.
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
# NOTE: the boundary(...) term needs the patched Finch (boundary integrals in
# FEM weak forms did not compile on the released code, see notes).
#
# Continuous Galerkin, linear elements on a quad mesh.
#==============================================================================#

using Finch

initFinch("cdi2d");

# results are written to output/, a folder kept in the repo
useLog("cdi2dlog", dir="output", level=3)

domain(2)
functionSpace(order=1)
timeStepper(EULER_EXPLICIT)

# The system matrix is the (symmetric positive definite) mass matrix, so
# plain conjugate gradients converge in a few iterations from the previous
# solution. The direct solver is far slower for this many time steps.
linAlgOptions(iterative=true, method="CG", pc="NONE")

nel = 300
mesh(QUADMESH, elsperdim=nel)      # unit square

# Name the four walls. The built-in quad mesh gives the whole boundary ID 1
# (its bids option is ignored), so the sides are assigned here by position.
addBoundaryID(1, "x==0")           # left,   inflow
addBoundaryID(2, "x==1")           # right,  outflow
addBoundaryID(3, "y==0")           # bottom, u.n = 0
addBoundaryID(4, "y==1")           # top,    u.n = 0

phi  = variable("phi")     # unknown at the new time level
phio = variable("phio")    # known field: previous time level
nhat = variable("nhat", type=VECTOR)   # known field: nodal unit normal, built from phio
testSymbol("v")

# Physical parameters
dx   = 1.0/nel
ux   = 1.0                 # advection velocity components
uy   = 0.0
umax = sqrt(ux^2 + uy^2)
gam  = umax                # sharpening parameter gamma (= max |u|)
eps  = dx                  # interface thickness (eps/dx >= (u/gamma + 1)/2)
R    = 0.15                # circle radius
xc0  = 0.3                 # initial circle center
yc0  = 0.5

# Inflow: the fluid entering at x = 0 carries phi = 0, so that value is imposed.
# The other three walls get no strong condition (NO_BC, the default, so
# nothing is declared). What crosses them is set by the boundary integral in
# the weak form below; on the Dirichlet nodes that integral is discarded with
# the rest of the row, so it only acts at the outflow.
boundary(phi, 1, DIRICHLET, 0)     # left, inflow

# Initial condition: equilibrium tanh profile around the circle
ic = "0.5*(1 - tanh((sqrt((x-$xc0)^2 + (y-$yc0)^2) - $R)/(2*$eps)))"
initial(phi, ic)

# Time stepping. Explicit diffusive limit in 2D (lumped mass):
# dt <= dx^2/(4*gamma*eps). Consistent mass matrix needs roughly 1/3 of that.
T  = 0.4
dt = 0.2*dx*dx/(4*gam*eps)
nsteps = Int(round(T/dt))
setSteps(dt, nsteps)

# Coefficients (theta as in the derivation, dimensional form)
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", [string(ux), string(uy)], type=VECTOR)   # velocity vector
coefficient("gamma", gam)
coefficient("epsilon", eps)
# (the 1e-12 that regularizes |grad(phi)| in the bulk is passed to recover_nhat! below)

# Nodal interface normal by finite element gradient recovery (lumped L2
# projection). The element gradient of the linear solution is discontinuous
# between elements; its projection onto the nodal basis,
#     sum_e int_e g_i psi_l = sum_e int_e (dphi_h/dx_i) psi_l,
# with the mass matrix lumped, is the area-weighted average of the element
# gradients around each node and is continuous. It uses only Finch's element
# connectivity, reference element and Jacobians, so it works on any mesh; on
# this uniform grid it reduces to a central difference smoothed over the
# neighbouring rows. The normal is then formed at the nodes.
const g_grid  = Finch.finch_state.grid_data
const g_refel = Finch.finch_state.refel
const g_geo   = Finch.finch_state.geo_factors
const g_Np    = g_refel.Np
const g_Nqp   = g_refel.Nqp
const g_nnode = size(g_grid.allnodes, 2)
const g_nel   = size(g_grid.loc2glb, 2)
# lumped mass: int psi_l, assembled once
const g_lump = zeros(g_nnode)
for e = 1:g_nel
    detj = g_geo.detJ[1, e]
    for i = 1:g_Np
        acc = 0.0
        for q = 1:g_Nqp
            acc += g_refel.wg[q] * g_refel.Q[q, i]
        end
        g_lump[g_grid.loc2glb[i, e]] += acc * detj
    end
end
const g_gx = zeros(g_nnode)
const g_gy = zeros(g_nnode)
const g_pe = zeros(g_Np)
const g_dr = zeros(g_Nqp)
const g_ds = zeros(g_Nqp)

# Project the gradient of the nodal field p onto the nodes, then normalize.
function recover_nhat!(nhatvals, p, delta)
    fill!(g_gx, 0.0); fill!(g_gy, 0.0)
    @inbounds for e = 1:g_nel
        J = g_geo.J[e]                      # constant Jacobian on this mesh: length-1 vectors
        rx = J.rx[1]; ry = J.ry[1]; sx = J.sx[1]; sy = J.sy[1]
        detj = g_geo.detJ[1, e]
        for i = 1:g_Np
            g_pe[i] = p[g_grid.loc2glb[i, e]]
        end
        # dphi/dr, dphi/ds at the quadrature points
        for q = 1:g_Nqp
            dr = 0.0; ds = 0.0
            for i = 1:g_Np
                dr += g_refel.Qr[q, i] * g_pe[i]
                ds += g_refel.Qs[q, i] * g_pe[i]
            end
            g_dr[q] = dr; g_ds[q] = ds
        end
        # int (dphi/dx) psi_l and int (dphi/dy) psi_l over the element
        for i = 1:g_Np
            ax = 0.0; ay = 0.0
            for q = 1:g_Nqp
                w  = g_refel.wg[q] * g_refel.Q[q, i]
                ax += w * (g_dr[q]*rx + g_ds[q]*sx)
                ay += w * (g_dr[q]*ry + g_ds[q]*sy)
            end
            k = g_grid.loc2glb[i, e]
            g_gx[k] += ax * detj
            g_gy[k] += ay * detj
        end
    end
    @inbounds for k = 1:g_nnode
        gx = g_gx[k] / g_lump[k]
        gy = g_gy[k] / g_lump[k]
        m  = sqrt(gx*gx + gy*gy + delta)
        nhatvals[1, k] = gx/m
        nhatvals[2, k] = gy/m
    end
    return nothing
end

# g is the normal component of the total flux on the boundary, as in the
# derivation:  g = (theta2*u*phi - theta3*gamma*eps*grad(phi)
#                   + theta4*gamma*(phi - phi^2)*nhat) . n
# At an outflow boundary only the convective part is kept and it is evaluated
# from the previous level rather than prescribed; the diffusive and
# sharpening parts are taken as zero. normal() is the outward unit normal.
parameter("g", "theta2*phio*dot(u, normal())")

# Snapshots for the time series. Saving all 2400 steps would be wasteful, so
# phi is kept every 40th step: 60 snapshots, plus one at t = 0 further down,
# makes the 61 the MATLAB scripts read.
snap_every = 40

# A list of lists. Each entry is a copy of every nodal value, so snapshots[k] is
# the whole field at snapshot k.
snapshots  = Vector{Vector{Float64}}()

# The simulation time of each snapshot.
snap_times = Float64[]

# Step counter. A Ref is a box holding the value. The post-step body is its own
# function, so a plain Int assigned in there would just become a local variable;
# writing step_count[] updates this one.
step_count = Ref(0)

# After each step the new solution becomes the "previous" field, and a
# snapshot is stored when due.
# Must be declared before weakForm (code is generated there).
@postStepFunction(
    begin
        phio.values .= phi.values
        recover_nhat!(nhat.values, phio.values, 1e-12)   # nodal normal from the new previous level
        step_count[] += 1
        if step_count[] % snap_every == 0
            push!(snapshots, copy(phi.values[:]))
            push!(snap_times, step_count[] * dt)
        end
    end
)

# Weak form (residual = 0), all spatial terms at the previous level phio:
#   theta1*(Dt phi, v)
#   - theta2/2*(u*phio, grad v)                       conservative half
#   - theta2/2*(phio, div(u)*v + u.grad v)            non-conservative half
#   + theta3*(gamma*eps*grad phio, grad v)            diffusion
#   - theta4*(gamma*(phio - phio^2)*nhat, grad v)     sharpening, nodal nhat
#   + <g, v>_boundary                                 convective outflow
# with nhat and g the parameters defined above. The boundary term is what is
# left over from integrating the two advection halves by parts (each gives
# half of it); the diffusive and sharpening boundary terms are dropped. It
# uses phio, so it only adds to the right-hand side.
weakForm(phi, "theta1*Dt(phi)*v" *
              " - 0.5*theta2*phio*dot(u, grad(v))" *
              " - 0.5*theta2*phio*(div(u)*v + dot(u, grad(v)))" *
              " + theta3*gamma*epsilon*dot(grad(phio), grad(v))" *
              " - theta4*gamma*(phio - phio*phio) .* dot(nhat, grad(v))" *
              " + boundary(g*v)")

# Evaluate the initial condition now so it can be stored as the t = 0 snapshot
evalInitialConditions()
phio.values .= phi.values   # phio holds the field the first step reads
recover_nhat!(nhat.values, phio.values, 1e-12)
push!(snapshots, copy(phi.values[:]))
push!(snap_times, 0.0)
# Lumped finite element weights: the lumped mass entry of a node is the integral
# of its basis function, which the boundary truncates to half at an edge node and
# a quarter at a corner. This keeps the mass sum exact once phi reaches the wall.
nodes = Finch.finch_state.grid_data.allnodes
wts = fill(dx*dx, size(nodes, 2))
for i = 1:size(nodes, 2)
    if nodes[1,i] < 1e-12 || nodes[1,i] > 1 - 1e-12; wts[i] *= 0.5; end
    if nodes[2,i] < 1e-12 || nodes[2,i] > 1 - 1e-12; wts[i] *= 0.5; end
end
mass0 = sum(wts .* phi.values[:])

solve(phi)

finalizeFinch()

# Error check against the translated exact profile
x = Finch.finch_state.grid_data.allnodes[1, :]
y = Finch.finch_state.grid_data.allnodes[2, :]
phivals = phi.values[:]
xc = xc0 + ux*T
yc = yc0 + uy*T
r  = sqrt.((x .- xc).^2 .+ (y .- yc).^2)
exact = 0.5 .* (1 .- tanh.((r .- R) ./ (2*eps)))

maxerr = maximum(abs.(phivals .- exact))
massT  = sum(wts .* phivals)
println("max error = " * string(maxerr))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))
println("relative mass change = " * string((massT - mass0)/mass0))

using DelimitedFiles

data = [x y phivals exact]
writedlm("output/solution.csv", data, ',')

# Time series: columns 1-2 = x, y; columns 3.. = phi at each snapshot time
ts = hcat(x, y, snapshots...)
writedlm("output/solution_timeseries.csv", ts, ',')
writedlm("output/snapshot_times.csv", snap_times, ',')
println(string(length(snap_times)) * " snapshots written")
