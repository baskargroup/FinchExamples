#==============================================================================#
# 1D Conservative Diffusive Interface (CDI) model, solved explicitly with
# Finch V2 (branch dev). Run with ./run.sh; the target (CPU Julia, or MLIR on
# the CPU/GPU) is chosen by FINCH_TARGET, see ../common.jl.
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

include(joinpath(@__DIR__, "..", "common.jl"))   # Finch, FINCH_TARGET switch, helpers

# results are written to output/, a folder kept in the repo; Finch's generated
# project (and its raw results) goes to output/finch/
const OUT_DIR  = joinpath(@__DIR__, get(ENV, "CDI_OUTDIR", "output"))   # CDI_OUTDIR: another results folder
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
const NSTEPS = Int(round(T/(0.2*DX*DX/(2*GAM*EPS))))
const DT     = T/NSTEPS           # the step Finch takes (Nsteps is what it is given)

# Mass matrices only (the projection and the explicit step): CG
ctx = init_finch("cdi1d"; finch_solver_kwargs(:cg)...,
                 target = finch_target("cdi1d", PROJ_DIR))

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
time_stepper(STEP_EULER_EXPLICIT; Nsteps = NSTEPS, store_trajectory = true, snapshot_every = 25)

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
               " - theta4*gamma*(phi - phi*phi)*" * NDOTGV_NODAL_1D * "*grad(v)")

println("dt = $DT, steps = $NSTEPS")
step_schedule(gx, phi)
solve()
finch_run!(ctx, phi)

# Nodes in order of x
x = ctx.grid.nodes[1, :]
o = sortperm(x)
x = x[o]

# Compare with the exact profile at the final time (the last snapshot: both targets write it)
snap_times, fields = read_snapshots(PROJ_DIR)
phivals = final_values(fields, "phi")[o]
exact = 0.5 .* (1 .+ tanh.((x .- X0 .- UVEL*T) ./ (2*EPS)))

maxerr = maximum(abs.(phivals .- exact))
println("max error = " * string(maxerr))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))

data = [x phivals exact]
writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')

# Time series: column 1 is x, the rest is phi at each snapshot time
write_timeseries(OUT_DIR, "solution", reshape(x, :, 1), fields["phi"][o, :])
write_snapshot_times(OUT_DIR, snap_times)
println(string(length(snap_times)) * " snapshots written")
