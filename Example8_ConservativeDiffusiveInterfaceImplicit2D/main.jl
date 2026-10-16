#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, semi-implicit, in Finch.
#
#   theta1*phi_t + theta2*(1/2)*[div(u*phi) + u.grad(phi)]
#     - theta3*div(gamma*eps*grad(phi)) + theta4*div(gamma*(phi - phi^2)*nhat) = 0
#
# phi is the phase field, 0 in one phase and 1 in the other, and nhat is the
# unit interface normal grad(phi)/|grad(phi)|. The diffusion and sharpening
# terms act against each other and hold the interface at a tanh profile of
# fixed width.
#
# Time discretization (semi-implicit, two Picard iterations per step):
#   - BDF2 for the time derivative, BDF1 on the first step:
#       phi_t ~ (a0*phi + a1*phio + a2*phioo)/dt
#   - advection and diffusion implicit (linear in phi)
#   - sharpening linearized with a lagged field phis:
#       gamma*(phi - phi^2)*nhat ~ gamma*phi*(1 - phis)*nhat(phis)
#   - iteration 1 uses phis = phio, iteration 2 uses phis = the result of
#     iteration 1. Each iteration is one linear solve; the matrix depends on
#     phis through the sharpening term, so it is rebuilt for iteration 2.
#
# Finch's built-in time loop cannot run two solves per step with the lagged
# field updated in between, so the time loop is written here in Julia and
# Finch is given a steady linear problem for one iteration (see Example 7).
# The matrix is not symmetric because of advection, so GMRES is used.
#
# Test problem: a circle of radius R carried by a uniform velocity, the same
# case as Example 6. The profile translates with almost no change of shape, so
#   phi(x,y,t) = 0.5*(1 - tanh((r - R)/(2*eps)))
# measured from the moving centre is close to the solution; a curved interface
# is not an equilibrium profile, so it is not exact.
#
# Boundaries: phi = 0 is imposed at the inflow x = 0. At the outflow x = 1 the
# convective boundary flux g = theta2*(u.n)*phi is kept in the weak form and the
# diffusive and sharpening fluxes are dropped (the do-nothing outflow
# condition), so mass leaves at the rate the flow carries it out. g is
# evaluated from the lagged field phis and enters the right-hand side, which is
# where the derivation places it; with phis updated in each Picard iteration it
# approaches the implicit value. Top and bottom have u.n = 0.
#
# Continuous Galerkin, linear elements on a quad mesh.
#==============================================================================#

using Finch

initFinch("cdi2dimplicit");

# results are written to output/, a folder kept in the repo
useLog("cdi2dimplicitlog", dir="output", level=2)

domain(2)
functionSpace(order=1)
# No timeStepper() call: Finch sees a steady problem. Time stepping is below.

# The system matrix changes every Picard iteration and is not symmetric
# (advection), so it is assembled each time and solved with GMRES. The default
# relative tolerance of 1e-8 shows up directly as a mass drift of that size;
# 1e-12 removes it at almost no extra cost (about 0.35 s per solve either way).
linAlgOptions(iterative=true, method="GMRES", pc="NONE", reltol=1e-12)

nel = 300
mesh(QUADMESH, elsperdim=nel)      # unit square

# Name the four walls. The built-in quad mesh gives the whole boundary ID 1
# (its bids option is ignored), so the sides are assigned here by position.
addBoundaryID(1, "x==0")           # left,   inflow
addBoundaryID(2, "x==1")           # right,  outflow
addBoundaryID(3, "y==0")           # bottom, u.n = 0
addBoundaryID(4, "y==1")           # top,    u.n = 0

phi   = variable("phi")    # unknown: the current Picard iterate
phio  = variable("phio")   # known: previous time level
phioo = variable("phioo")  # known: two time levels back
phis  = variable("phis")   # known: lagged field for sharpening and outflow
testSymbol("v")

# Physical parameters (same case as Example 6)
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

# Time step. Diffusion is implicit, so the explicit limit dx^2/(4*gamma*eps)
# does not apply; dt is set by accuracy through the advective CFL number.
# CFL 0.25 kept phi within [0,1] and sat near the spatial error floor in the
# 1D study (Example 7).
T   = 0.4
CFL = 0.25
dt  = CFL*dx/umax
nsteps = Int(round(T/dt))

# BDF coefficients, divided by dt. They are read by Finch at every assembly
# through the coefficients a0dt, a1dt, a2dt below, so changing these numbers
# between solves is enough to switch from BDF1 (first step) to BDF2.
a0_dt = 1.0/dt;  a1_dt = -1.0/dt;  a2_dt = 0.0           # BDF1
bdf2!() = (global a0_dt = 1.5/dt; global a1_dt = -2.0/dt; global a2_dt = 0.5/dt)

# Coefficients (theta as in the derivation, dimensional form). The a*dt
# coefficients are string expressions that evaluate the Julia globals above
# (Main. is needed because Finch evaluates them inside its own module).
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", [string(ux), string(uy)], type=VECTOR)   # velocity vector
coefficient("gamma", gam)
coefficient("epsilon", eps)
coefficient("delta", 1e-12)                     # keeps nhat finite where grad(phi) = 0
coefficient("a0dt", "Main.a0_dt")
coefficient("a1dt", "Main.a1_dt")
coefficient("a2dt", "Main.a2_dt")

# The lagged unit interface normal, built from phis
parameter("nhat", "grad(phis) ./ sqrt.(dot(grad(phis), grad(phis)) .+ delta)")

# g is the normal component of the total flux on the boundary, as in the
# derivation:  g = (theta2*u*phi - theta3*gamma*eps*grad(phi)
#                   + theta4*gamma*(phi - phi^2)*nhat) . n
# At an outflow boundary only the convective part is kept, evaluated from the
# lagged field; the diffusive and sharpening parts are taken as zero.
# normal() is the outward unit normal.
parameter("g", "theta2*phis*dot(u, normal())")

# Weak form of one Picard iteration, written as a residual that Finch sets to
# zero. Terms with phi go to the matrix, the rest to the right-hand side:
#
#   theta1*(a0/dt*phi, v) + theta1*((a1*phio + a2*phioo)/dt, v)   time derivative
#   - theta2/2*(u*phi, grad v)                                    advection, conservative half
#   - theta2/2*(phi, div(u)*v + u.grad v)                         advection, other half
#   + theta3*(gamma*eps*grad phi, grad v)                         diffusion
#   - theta4*(gamma*phi*(1 - phis)*nhat(phis), grad v)            sharpening, lagged
#   + <g, v>_boundary                                             convective outflow
weakForm(phi, "theta1*a0dt*phi*v" *
              " + theta1*(a1dt*phio + a2dt*phioo)*v" *
              " - 0.5*theta2*phi*dot(u, grad(v))" *
              " - 0.5*theta2*phi*(div(u)*v + dot(u, grad(v)))" *
              " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
              " - theta4*gamma*phi*(1 - phis) .* dot(nhat, grad(v))" *
              " + boundary(g*v)")

# Snapshots for the time series: phi every snap_every steps plus t = 0, which
# gives the same 61 columns the MATLAB scripts read.
snap_every = max(1, nsteps ÷ 60)
snapshots  = Vector{Vector{Float64}}()
snap_times = Float64[]

# Start: all history levels equal the initial condition
evalInitialConditions()
phio.values  .= phi.values
phioo.values .= phi.values
phis.values  .= phi.values
push!(snapshots, copy(phi.values[:]))
push!(snap_times, 0.0)

# Lumped finite element weights for the mass: the lumped mass entry of a node
# is the integral of its basis function, which the boundary truncates to half
# at an edge node and a quarter at a corner.
nodes = Finch.finch_state.grid_data.allnodes
wts = fill(dx*dx, size(nodes, 2))
for i = 1:size(nodes, 2)
    if nodes[1,i] < 1e-12 || nodes[1,i] > 1 - 1e-12; wts[i] *= 0.5; end
    if nodes[2,i] < 1e-12 || nodes[2,i] > 1 - 1e-12; wts[i] *= 0.5; end
end
mass0 = sum(wts .* phi.values[:])

# Time loop. Two Picard iterations per step, then shift the history.
t0 = time()
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
    if step % 10 == 0 || step == nsteps
        println("step $step / $nsteps   t = $(round(step*dt, digits=5))   elapsed $(round(time() - t0, digits=1)) s")
        flush(stdout)
    end
end

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
println("dt = " * string(dt) * ", steps = " * string(nsteps) * " (CFL " * string(CFL) * ")")
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
