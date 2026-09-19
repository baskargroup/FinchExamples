#==============================================================================#
# 1D Conservative Diffusive Interface (CDI) model, fully explicit, no VMS.
#
#   theta1*phi_t + theta2*(1/2)*[(u*phi)_x + u*phi_x]
#     - theta3*(gamma*eps*phi_x)_x + theta4*(gamma*(phi - phi^2)*nhat)_x = 0
#
# Fine scale phi' = 0 (no stabilization). Forward Euler in time: the unknown
# phi appears only in the mass term, everything else is evaluated at the
# previous time level, stored in the variable phio. One linear solve per step.
#
# Test: single interface advected by a uniform velocity. The equilibrium tanh
# profile translates without change of shape:
#   phi_exact(x,t) = 0.5*(1 + tanh((x - x0 - u*t)/(2*eps)))
# CG, linear elements.
#==============================================================================#

using Finch

initFinch("cdi1d");

useLog("cdi1dlog", level=3)

domain(1)
functionSpace(order=1)
timeStepper(EULER_EXPLICIT)

nel = 200
mesh(LINEMESH, elsperdim=nel, bids=2)

phi  = variable("phi")     # unknown at the new time level
phio = variable("phio")    # known field: previous time level
testSymbol("v")

# Physical parameters
dx    = 1.0/nel
uvel  = 1.0                # advection velocity
gam   = 1.0                # sharpening parameter gamma (= max |u|)
eps   = dx                 # interface thickness (eps/dx >= (u/gamma + 1)/2)
x0    = 0.25               # initial interface position

# Boundary conditions: phi = 0 on the left, phi = 1 on the right
boundary(phi, 1, DIRICHLET, 0)
boundary(phi, 2, DIRICHLET, 1)
boundary(phio, 1, DIRICHLET, 0)
boundary(phio, 2, DIRICHLET, 1)

# Initial condition: equilibrium tanh profile
ic = "0.5*(1 + tanh((x - $x0)/(2*$eps)))"
initial(phi, ic)
initial(phio, ic)

# Time stepping. Explicit limit (lumped mass): dt <= dx^2/(2*gamma*eps).
# Consistent mass matrix needs roughly 1/3 of that.
T  = 0.25
dt = 0.2*dx*dx/(2*gam*eps)
nsteps = Int(round(T/dt))
setSteps(dt, nsteps)

# Coefficients (theta as in the derivation, dimensional form)
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", [string(uvel)], type=VECTOR)   # velocity vector (1 component in 1D)
coefficient("gamma", gam)
coefficient("epsilon", eps)
coefficient("delta", 1e-12)                     # regularizes |grad(phi)| in the bulk

# Snapshots for the time series: store phi every snap_every steps
snap_every = 25
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
# Note: in 1D Finch treats a 1-component vector as a scalar and rejects div(u),
# so div(u) is written as grad(u) here (identical in 1D). In 2D use div(u).
#   + theta3*(gamma*eps*grad phio, grad v)            diffusion
#   - theta4*(gamma*(phio - phio^2)*nhat, grad v)     sharpening
# with nhat = grad(phio)/sqrt(|grad(phio)|^2 + delta)
weakForm(phi, "theta1*Dt(phi*v)" *
              " - 0.5*theta2*phio*dot(u, grad(v))" *
              " - 0.5*theta2*phio*(grad(u)*v + dot(u, grad(v)))" *
              " + theta3*gamma*epsilon*dot(grad(phio), grad(v))" *
              " - theta4*gamma*(phio - phio*phio) .* dot(grad(phio) ./ sqrt.(dot(grad(phio), grad(phio)) .+ delta), grad(v))")

# Evaluate the initial condition now so it can be stored as the t = 0 snapshot
evalInitialConditions()
push!(snapshots, copy(phi.values[:]))
push!(snap_times, 0.0)

solve(phi)

finalizeFinch()

# Error check against the translated exact profile
x = Finch.finch_state.grid_data.allnodes[1, :]
phivals = phi.values[:]
exact = 0.5 .* (1 .+ tanh.((x .- x0 .- uvel*T) ./ (2*eps)))

maxerr = maximum(abs.(phivals .- exact))
println("max error = " * string(maxerr))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))

using DelimitedFiles

data = [x phivals exact]
writedlm("solution.csv", data, ',')

# Time series: column 1 = x, columns 2.. = phi at each snapshot time
ts = hcat(x, snapshots...)
writedlm("solution_timeseries.csv", ts, ',')
writedlm("snapshot_times.csv", snap_times, ',')
println(string(length(snap_times)) * " snapshots written")
