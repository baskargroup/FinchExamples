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
coefficient("delta", 1e-12)                     # regularizes |grad(phi)| in the bulk

# The unit interface normal, named here so the weak form below stays readable.
# parameter() is not a coefficient: it is just a named symbolic expression that
# Finch substitutes into the weak form, so the generated code is unchanged.
parameter("nhat", "grad(phio) ./ sqrt.(dot(grad(phio), grad(phio)) .+ delta)")

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
#   - theta4*(gamma*(phio - phio^2)*nhat, grad v)     sharpening
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
