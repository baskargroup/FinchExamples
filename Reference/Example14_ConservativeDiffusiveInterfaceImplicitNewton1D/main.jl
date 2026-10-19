#==============================================================================#
# 1D Conservative Diffusive Interface (CDI) model, semi-implicit, in Finch,
# with the sharpening term linearized by Newton's method (compare Examples 8 and 10).
#
#   theta1*phi_t + theta2*(1/2)*[(u*phi)_x + u*phi_x]
#     - theta3*(gamma*eps*phi_x)_x + theta4*(gamma*(phi - phi^2)*nhat)_x = 0
#
# phi is the phase field, 0 in one phase and 1 in the other, and nhat is the
# unit interface normal phi_x/|phi_x|. The diffusion and sharpening terms act
# against each other and hold the interface at a tanh profile of fixed width.
#
# Time discretization (semi-implicit, two Picard iterations per step):
#   - BDF2 for the time derivative, BDF1 on the first step:
#       phi_t ~ (a0*phi + a1*phio + a2*phioo)/dt
#   - advection and diffusion implicit (linear in phi)
#   - sharpening linearized by Newton's method about a lagged field phis:
#       phi - phi^2 ~ phis*(1 - phis) + (1 - 2*phis)*(phi - phis)
#                   = (1 - 2*phis)*phi + phis^2
#       gamma*(phi - phi^2)*nhat ~ gamma*((1 - 2*phis)*phi + phis^2)*nhat(phis)
#     The part with phi goes to the matrix, phis^2 to the right-hand side; the
#     normal stays lagged. Example 8 keeps one factor implicit,
#     gamma*phi*(1 - phis)*nhat(phis); Example 10 lags the whole term. All three
#     agree once the iteration has converged (phi = phis); this one is the
#     consistent linearization, so the iteration converges quadratically.
#     The implicit coefficient gamma*(1 - 2*phis)*nhat is the sharpening
#     velocity: it vanishes at phi = 0.5 and changes sign across the interface.
#   - iteration 1 uses phis = phio, iteration 2 uses phis = the result of
#     iteration 1. Each iteration is one linear solve; the matrix depends on
#     phis, so matrix and right-hand side are both rebuilt.
#
# Finch's built-in time loop cannot run two solves per step with the lagged
# field updated in between, so the time loop is written here in Julia and
# Finch is given a steady linear problem for one iteration. Each call to
# solve(phi) reassembles the matrix and right-hand side from the current
# values of the known fields and coefficients.
#
# Test problem: one interface carried by a uniform velocity. The equilibrium
# profile translates without changing shape, so the exact solution is
#   phi(x,t) = 0.5*(1 + tanh((x - x0 - u*t)/(2*eps)))
#
# Continuous Galerkin, linear elements.
#==============================================================================#

using Finch

initFinch("cdi1dnewton");

# results are written to output/, a folder kept in the repo
useLog("cdi1dnewtonlog", dir="output", level=2)

domain(1)
functionSpace(order=1)
# No timeStepper() call: Finch sees a steady problem. Time stepping is below.

nel = 200
mesh(LINEMESH, elsperdim=nel, bids=2)

phi   = variable("phi")    # unknown: the current Picard iterate
phio  = variable("phio")   # known: previous time level
phioo = variable("phioo")  # known: two time levels back
phis  = variable("phis")   # known: lagged field for the sharpening term
testSymbol("v")

# Physical parameters (same case as Examples 5 and 8)
dx    = 1.0/nel
uvel  = 1.0                # advection velocity
gam   = 1.0                # sharpening parameter, usually max |u|
eps   = dx                 # interface thickness; eps/dx >= (u/gamma + 1)/2
                           # is what keeps phi within [0,1]
x0    = 0.25               # initial interface position

# The profile spans the domain, so the ends carry the two bulk values
boundary(phi, 1, DIRICHLET, 0)
boundary(phi, 2, DIRICHLET, 1)

# Initial condition: the equilibrium profile centred at x0
ic = "0.5*(1 + tanh((x - $x0)/(2*$eps)))"
initial(phi, ic)

# Time step. Diffusion is implicit, so the explicit limit dx^2/(2*gamma*eps)
# does not apply; dt is set by accuracy through the advective CFL number.
# CFL 0.25 keeps phi within [0,1] to round-off and sits near the spatial error
# floor; CFL 0.5 still beats the explicit scheme but overshoots by 4e-4, and
# CFL 1 and above lose accuracy quickly (see notes/).
T   = 0.25
CFL = 0.25
dt  = CFL*dx/uvel
nsteps = Int(round(T/dt))

# BDF coefficients, divided by dt. They are read by Finch at every assembly
# through the coefficients a0dt, a1dt, a2dt below, so changing these numbers
# between solves is enough to switch from BDF1 (first step) to BDF2.
a0_dt = 1.0/dt;  a1_dt = -1.0/dt;  a2_dt = 0.0           # BDF1
bdf2!() = (global a0_dt = 1.5/dt; global a1_dt = -2.0/dt; global a2_dt = 0.5/dt)

# Coefficients. The theta are the model coefficients of the derivation, all
# one in dimensional form. The a*dt coefficients are string expressions that
# evaluate the Julia globals above (Main. is needed because Finch evaluates
# them inside its own module).
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", [string(uvel)], type=VECTOR)   # velocity vector (1 component in 1D)
coefficient("gamma", gam)
coefficient("epsilon", eps)
coefficient("delta", 1e-12)                     # keeps nhat finite where grad(phi) = 0
coefficient("a0dt", "Main.a0_dt")
coefficient("a1dt", "Main.a1_dt")
coefficient("a2dt", "Main.a2_dt")

# The lagged unit interface normal, built from phis
parameter("nhat", "grad(phis) ./ sqrt.(dot(grad(phis), grad(phis)) .+ delta)")

# Weak form of one Picard iteration, written as a residual that Finch sets to
# zero. Terms with phi go to the matrix, the rest to the right-hand side:
#
#   theta1*(a0/dt*phi, v) + theta1*((a1*phio + a2*phioo)/dt, v)   time derivative
#   - theta2/2*(u*phi, grad v)                                    advection, conservative half
#   - theta2/2*(phi, div(u)*v + u.grad v)                         advection, other half
#   + theta3*(gamma*eps*grad phi, grad v)                         diffusion
#   - theta4*(gamma*((1 - 2*phis)*phi + phis^2)*nhat(phis), grad v)  sharpening, Newton
#
# The boundary terms from integrating by parts drop out because v vanishes at
# the two Dirichlet ends. In 1D Finch reads a one-component vector as a scalar
# and rejects div(u), so grad(u) appears instead; use div(u) in 2D.
weakForm(phi, "theta1*a0dt*phi*v" *
              " + theta1*(a1dt*phio + a2dt*phioo)*v" *
              " - 0.5*theta2*phi*dot(u, grad(v))" *
              " - 0.5*theta2*phi*(grad(u)*v + dot(u, grad(v)))" *
              " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
              " - theta4*gamma*((1 - 2*phis)*phi + phis*phis) .* dot(nhat, grad(v))")

# Snapshots for the time series: phi every snap_every steps plus t = 0, which
# gives the same 21 columns the MATLAB scripts read.
snap_every = max(1, nsteps ÷ 20)
snapshots  = Vector{Vector{Float64}}()
snap_times = Float64[]

# Start: all history levels equal the initial condition
evalInitialConditions()
phio.values  .= phi.values
phioo.values .= phi.values
phis.values  .= phi.values
push!(snapshots, copy(phi.values[:]))
push!(snap_times, 0.0)

# Time loop. Two Picard iterations per step, then shift the history.
for step = 1:nsteps
    if step == 2
        bdf2!()                          # BDF1 only on the first step
    end
    phis.values .= phio.values           # iteration 1: lag with the previous level
    solve(phi)
    phis.values .= phi.values            # iteration 2: lag with iteration 1
    solve(phi)
    phioo.values .= phio.values          # shift history
    phio.values  .= phi.values
    if step % snap_every == 0
        push!(snapshots, copy(phi.values[:]))
        push!(snap_times, step*dt)
    end
end

finalizeFinch()

# Compare with the exact profile at the final time
x = Finch.finch_state.grid_data.allnodes[1, :]
phivals = phi.values[:]
exact = 0.5 .* (1 .+ tanh.((x .- x0 .- uvel*T) ./ (2*eps)))

maxerr = maximum(abs.(phivals .- exact))
println("dt = " * string(dt) * ", steps = " * string(nsteps) * " (CFL " * string(CFL) * ")")
println("max error = " * string(maxerr))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))

using DelimitedFiles

data = [x phivals exact]
writedlm("output/solution.csv", data, ',')

# Time series: column 1 is x, the rest is phi at each snapshot time
ts = hcat(x, snapshots...)
writedlm("output/solution_timeseries.csv", ts, ',')
writedlm("output/snapshot_times.csv", snap_times, ',')
println(string(length(snap_times)) * " snapshots written")
