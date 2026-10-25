#==============================================================================#
# 1D Conservative Diffusive Interface (CDI) model, solved explicitly with
# Finch V2 (Julia target).
#
#   theta1*phi_t + theta2*(1/2)*[(u*phi)_x + u*phi_x]
#     - theta3*(gamma*eps*phi_x)_x + theta4*(gamma*(phi - phi^2)*nhat)_x = 0
#
# phi is the phase field, 0 in one phase and 1 in the other, and nhat is the
# unit interface normal phi_x/|phi_x|. The diffusion and sharpening terms act
# against each other and hold the interface at a tanh profile of fixed width.
#
# Interface normal: recovered at the NODES by an L2 projection of the finite
# element gradient of the previous level, (g, v) = (phi_x, v), a field solved
# every step before phi; nhat = g/sqrt(g^2 + delta). Finch V2 has no lumped
# mass matrix, so the projection uses the consistent one (the V1 version of
# this example lumped it). Taking the normal element by element from the
# gradient of the linear basis gives a field that is constant per element and
# jumps between elements; in 2D that produced grid-scale wiggles on interfaces
# running along the grid (see Examples 6 and 7). In 1D the two agree except in
# the tails, and the nodal form is used here for consistency.
#
# Forward Euler in time: every spatial term is evaluated at the previous level
# (Finch V2 does this itself under STEP_EULER_EXPLICIT), so each step is one
# solve with the constant mass matrix.
#
# Test problem: one interface carried by a uniform velocity. The equilibrium
# profile translates without changing shape, so the exact solution is
#   phi(x,t) = 0.5*(1 + tanh((x - x0 - u*t)/(2*eps)))
#
# Continuous Galerkin, linear elements.
#==============================================================================#

using Finch
using DelimitedFiles

# results are written to output/, a folder kept in the repo; Finch's generated
# project (and its raw results) goes to output/finch/
const OUT_DIR  = joinpath(@__DIR__, "output")
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

# Physical parameters
const NEL  = 200
const DX   = 1.0/NEL
const UVEL = 1.0                # advection velocity
const GAM  = 1.0                # sharpening parameter, usually max |u|
const EPS  = DX                 # interface thickness; eps/dx >= (u/gamma + 1)/2
                                # is what keeps phi within [0,1]
const X0   = 0.25               # initial interface position

# Time step. The explicit diffusive limit is dt <= dx^2/(2*gamma*eps) with a
# lumped mass matrix; Finch uses the consistent one, which needs about a third
# of that, so the factor below is well inside the stable range.
const T      = 0.25
const DT     = 0.2*DX*DX/(2*GAM*EPS)
const NSTEPS = Int(round(T/DT))

# Mass matrices only (the projection and the explicit step): CG
ctx = init_finch("cdi1d"; solver = :cg,
                 target = JuliaCPUTarget(out_dir = PROJ_DIR, func_name = "cdi1d"))

domain(1)
mesh(MESH_LINE; elsperdim = NEL)
function_space(order = 1)

phi = variable("phi")       # the phase field
gx  = variable("gx")        # projected gradient of phi, for the nodal normal
test_symbol("v")

# The profile spans the domain, so the ends carry the two bulk values
boundary(phi, 1, BC_DIRICHLET, 0.0)
boundary(phi, 2, BC_DIRICHLET, 1.0)

# Initial condition: the equilibrium profile centred at x0
initial(phi, (x, y, z, t) -> 0.5*(1 + tanh((x - X0)/(2*EPS))))

# Coefficients. The theta are the model coefficients of the derivation, all
# one in dimensional form.
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", UVEL)        # velocity (one component in 1D)
coefficient("gamma", GAM)
coefficient("epsilon", EPS)
coefficient("delta", 1e-12)   # keeps nhat finite where grad(phi) = 0

# Snapshots for the time series: phi every 25th step, 20 snapshots plus the
# one at t = 0, the 21 the MATLAB scripts read
time_interval(T)
time_stepper(STEP_EULER_EXPLICIT; dt = DT, store_trajectory = true, snapshot_every = 25)

# Weak forms, written as residuals that Finch sets to zero. The projection:
#   (g, v) - (grad phi, v)
# and the phase field, every spatial term at the previous level:
#   theta1*(Dt phi, v)
#   - theta2/2*(u*phi, grad v)                       advection, conservative half
#   + theta2/2*(u*grad phi, v)                       advection, other half
#   + theta3*(gamma*eps*grad phi, grad v)            diffusion
#   - theta4*(gamma*(phi - phi^2)*nhat, grad v)      sharpening
# The two advection halves together are the skew-symmetric form. Integration by
# parts moved a derivative onto v in each term, and the boundary terms drop out
# because v vanishes at the two Dirichlet ends.
weak_form(gx, "gx*v - grad(phi)*v")
weak_form(phi, "theta1*Dt(phi)*v" *
               " - 0.5*theta2*u*phi*grad(v)" *
               " + 0.5*theta2*u*grad(phi)*v" *
               " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
               " - theta4*gamma*(phi - phi*phi)*gx*pow(gx*gx + delta, -0.5)*grad(v)")

step_schedule(gx, phi)
solve()
Finch.run_emitted_julia!(ctx.target, phi)

# Nodes in order of x
x = ctx.grid.nodes[1, :]
o = sortperm(x)
x = x[o]

# Compare with the exact profile at the final time
phivals = phi.values[1, o]
exact = 0.5 .* (1 .+ tanh.((x .- X0 .- UVEL*T) ./ (2*EPS)))

maxerr = maximum(abs.(phivals .- exact))
println("max error = " * string(maxerr))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))

data = [x phivals exact]
writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')

# Time series: column 1 is x, the rest is phi at each snapshot time
results = joinpath(PROJ_DIR, "results")
_, snaps = Finch.read_results_manifest(results)
snapshots  = [Finch.load_snapshot_values(results, s)["phi"][o] for s in snaps]
snap_times = [Float64(s["t"]) for s in snaps]
writedlm(joinpath(OUT_DIR, "solution_timeseries.csv"), hcat(x, snapshots...), ',')
writedlm(joinpath(OUT_DIR, "snapshot_times.csv"), snap_times, ',')
println(string(length(snap_times)) * " snapshots written")
