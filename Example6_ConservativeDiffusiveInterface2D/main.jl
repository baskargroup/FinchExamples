#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, fully explicit, no VMS.
#
#   theta1*phi_t + theta2*(1/2)*[div(u*phi) + u.grad(phi)]
#     - theta3*div(gamma*eps*grad(phi)) + theta4*div(gamma*(phi - phi^2)*nhat) = 0
#
# Fine scale phi' = 0 (no stabilization). Forward Euler in time: the unknown
# phi appears only in the mass term, everything else is evaluated at the
# previous time level, stored in the variable phio. One linear solve per step.
#
# Test: circle of radius R advected by a uniform velocity (ux, uy). The
# equilibrium profile translates without change of shape:
#   phi_exact(x,y,t) = 0.5*(1 - tanh((r - R)/(2*eps))),  r = |x - xc(t)|
# CG, linear elements on a quad mesh. No-flux (natural) boundaries: the weak
# form omits the boundary integral, which is the g = 0 condition, so mass is
# conserved to round-off. Dirichlet phi = 0 would instead discard the little
# mass that reaches the wall once the interface tail gets there.
#==============================================================================#

using Finch

initFinch("cdi2d");

useLog("cdi2dlog", level=3)

domain(2)
functionSpace(order=1)
timeStepper(EULER_EXPLICIT)

# The system matrix is the (symmetric positive definite) mass matrix, so
# plain conjugate gradients converge in a few iterations from the previous
# solution. The direct solver is far slower for this many time steps.
linAlgOptions(iterative=true, method="CG", pc="NONE")

nel = 300
mesh(QUADMESH, elsperdim=nel)      # unit square, single boundary ID

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

# Boundary conditions: no flux on the whole boundary (natural condition)
boundary(phi, 1, NO_BC)
boundary(phio, 1, NO_BC)

# Initial condition: equilibrium tanh profile around the circle
ic = "0.5*(1 - tanh((sqrt((x-$xc0)^2 + (y-$yc0)^2) - $R)/(2*$eps)))"
initial(phi, ic)
initial(phio, ic)

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

# Snapshots for the time series: store phi every snap_every steps
snap_every = 40
snapshots  = Vector{Vector{Float64}}()   # one vector of nodal values per snapshot
snap_times = Float64[]
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
# with nhat = grad(phio)/sqrt(|grad(phio)|^2 + delta)
weakForm(phi, "theta1*Dt(phi*v)" *
              " - 0.5*theta2*phio*dot(u, grad(v))" *
              " - 0.5*theta2*phio*(div(u)*v + dot(u, grad(v)))" *
              " + theta3*gamma*epsilon*dot(grad(phio), grad(v))" *
              " - theta4*gamma*(phio - phio*phio) .* dot(grad(phio) ./ sqrt.(dot(grad(phio), grad(phio)) .+ delta), grad(v))")

# Evaluate the initial condition now so it can be stored as the t = 0 snapshot
evalInitialConditions()
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
writedlm("solution.csv", data, ',')

# Time series: columns 1-2 = x, y; columns 3.. = phi at each snapshot time
ts = hcat(x, y, snapshots...)
writedlm("solution_timeseries.csv", ts, ',')
writedlm("snapshot_times.csv", snap_times, ',')
println(string(length(snap_times)) * " snapshots written")
