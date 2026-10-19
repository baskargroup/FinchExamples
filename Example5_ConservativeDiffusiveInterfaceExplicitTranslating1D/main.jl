#==============================================================================#
# 1D Conservative Diffusive Interface (CDI) model, solved explicitly in Finch.
#
#   theta1*phi_t + theta2*(1/2)*[(u*phi)_x + u*phi_x]
#     - theta3*(gamma*eps*phi_x)_x + theta4*(gamma*(phi - phi^2)*nhat)_x = 0
#
# phi is the phase field, 0 in one phase and 1 in the other, and nhat is the
# unit interface normal phi_x/|phi_x|. The diffusion and sharpening terms act
# against each other and hold the interface at a tanh profile of fixed width.
#
# Interface normal: recovered at the NODES by a lumped L2 projection of the
# finite element gradient of the previous level (the area-weighted average of
# the element gradients around each node), stored in a known variable and
# refreshed after every step. Taking the normal element by element from the
# gradient of the linear basis gives a field that is constant per element and
# jumps between elements; in 2D that produced grid-scale wiggles on interfaces
# running along the grid (see Examples 6 and 7). In 1D the two agree except in
# the tails, and the nodal form is used here for consistency.
#
# Forward Euler in time. The unknown phi then appears only in the time term,
# so each step is one linear solve with the constant mass matrix.
# Every spatial term is evaluated at the previous level, held in phio.
#
# Test problem: one interface carried by a uniform velocity. The equilibrium
# profile translates without changing shape, so the exact solution is
#   phi(x,t) = 0.5*(1 + tanh((x - x0 - u*t)/(2*eps)))
#
# Continuous Galerkin, linear elements.
#==============================================================================#

using Finch

initFinch("cdi1d");

# results are written to output/, a folder kept in the repo
useLog("cdi1dlog", dir="output", level=3)

domain(1)
functionSpace(order=1)
timeStepper(EULER_EXPLICIT)

nel = 200
mesh(LINEMESH, elsperdim=nel, bids=2)

phi  = variable("phi")     # unknown at the new time level
phio = variable("phio")    # known field: previous time level
nhat = variable("nhat")    # known field: nodal unit normal (one component in 1D), built from phio
testSymbol("v")

# Physical parameters
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
# Stores the IC expression, does not evaluate it yet.
ic = "0.5*(1 + tanh((x - $x0)/(2*$eps)))"
initial(phi, ic)

# Time step. The explicit diffusive limit is dt <= dx^2/(2*gamma*eps) with a
# lumped mass matrix; Finch uses the consistent one, which needs about a third
# of that, so the factor below is well inside the stable range.
T  = 0.25
dt = 0.2*dx*dx/(2*gam*eps)
nsteps = Int(round(T/dt))
setSteps(dt, nsteps)

# Coefficients. The theta are the model coefficients of the derivation, all
# one in dimensional form.
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", [string(uvel)], type=VECTOR)   # velocity vector (1 component in 1D)
coefficient("gamma", gam)
coefficient("epsilon", eps)
# (the 1e-12 that keeps nhat finite where grad(phi) = 0 is passed to recover_nhat! below)

# Nodal interface normal by finite element gradient recovery (lumped L2
# projection): for every node l,  sum_e int_e g psi_l = sum_e int_e (dphi_h/dx) psi_l
# with the mass matrix lumped. Uses Finch's element connectivity, reference
# element and Jacobians, so it is the same on any mesh.
const g_grid  = Finch.finch_state.grid_data
const g_refel = Finch.finch_state.refel
const g_geo   = Finch.finch_state.geo_factors
const g_Np    = g_refel.Np
const g_Nqp   = g_refel.Nqp
const g_nnode = size(g_grid.allnodes, 2)
const g_nel   = size(g_grid.loc2glb, 2)
const g_lump  = zeros(g_nnode)            # lumped mass: int psi_l, assembled once
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
const g_pe = zeros(g_Np)
const g_dr = zeros(g_Nqp)

# Project dphi/dx of the nodal field p onto the nodes, then normalize.
function recover_nhat!(nhatvals, p, delta)
    fill!(g_gx, 0.0)
    @inbounds for e = 1:g_nel
        rx = g_geo.J[e].rx[1]                 # constant Jacobian on this mesh
        detj = g_geo.detJ[1, e]
        for i = 1:g_Np
            g_pe[i] = p[g_grid.loc2glb[i, e]]
        end
        for q = 1:g_Nqp                       # dphi/dr at the quadrature points
            dr = 0.0
            for i = 1:g_Np
                dr += g_refel.Qr[q, i] * g_pe[i]
            end
            g_dr[q] = dr
        end
        for i = 1:g_Np                        # int (dphi/dx) psi_l over the element
            ax = 0.0
            for q = 1:g_Nqp
                ax += g_refel.wg[q] * g_refel.Q[q, i] * g_dr[q] * rx
            end
            g_gx[g_grid.loc2glb[i, e]] += ax * detj
        end
    end
    @inbounds for k = 1:g_nnode
        gx = g_gx[k] / g_lump[k]
        nhatvals[k] = gx / sqrt(gx*gx + delta)
    end
    return nothing
end

# Snapshots for the time series. Saving all 500 steps would be wasteful, so phi
# is kept every 25th step: 20 snapshots, plus one at t = 0 further down, makes
# the 21 the MATLAB scripts read.
snap_every = 25

# A list of lists. Each entry is a copy of every nodal value, so snapshots[k] is
# the whole field at snapshot k.
snapshots  = Vector{Vector{Float64}}()

# The simulation time of each snapshot.
snap_times = Float64[]

# Step counter. A Ref is a box holding the value. The post-step body is its own
# function, so a plain Int assigned in there would just become a local variable;
# writing step_count[] updates this one.
step_count = Ref(0)

# After each step the new solution becomes the previous field, and a snapshot
# is kept when due. This must come before weakForm, because that call is where
# Finch generates the time loop.
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

# Weak form, written as a residual that Finch sets to zero. Only the time term
# carries the unknown; every spatial term uses the previous level phio:
#
#   theta1*(Dt phi, v)
#   - theta2/2*(u*phio, grad v)                      advection, conservative half
#   - theta2/2*(phio, div(u)*v + u.grad v)           advection, other half
#   + theta3*(gamma*eps*grad phio, grad v)           diffusion
#   - theta4*(gamma*(phio - phio^2)*nhat, grad v)    sharpening
#
# The two advection halves together are the skew-symmetric form. Integration by
# parts moved a derivative onto v in each term, and the boundary terms drop out
# because v vanishes at the two Dirichlet ends.
#
# nhat is the parameter defined above. One thing to note in the code below: in
# 1D Finch reads a one-component vector as a scalar and rejects div(u), so
# grad(u) appears instead. The two are the same thing in 1D, but use div(u)
# in 2D.
weakForm(phi, "theta1*Dt(phi)*v" *
              " - 0.5*theta2*phio*dot(u, grad(v))" *
              " - 0.5*theta2*phio*(grad(u)*v + dot(u, grad(v)))" *
              " + theta3*gamma*epsilon*dot(grad(phio), grad(v))" *
              " - theta4*gamma*(phio - phio*phio) .* nhat*grad(v)")

# Evaluate the initial condition now so it can be kept as the t = 0 snapshot
evalInitialConditions()
phio.values .= phi.values   # phio holds the field the first step reads
recover_nhat!(nhat.values, phio.values, 1e-12)
push!(snapshots, copy(phi.values[:]))
push!(snap_times, 0.0)

solve(phi)

finalizeFinch()

# Compare with the exact profile at the final time
x = Finch.finch_state.grid_data.allnodes[1, :]
phivals = phi.values[:]
exact = 0.5 .* (1 .+ tanh.((x .- x0 .- uvel*T) ./ (2*eps)))

maxerr = maximum(abs.(phivals .- exact))
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
