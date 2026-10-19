#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, solved explicitly in Finch
# with VMS stabilization: Zalesak's slotted disk in solid-body rotation
# (compare Example 7, the same case without stabilization).
#
#   theta1*phi_t + theta2*(1/2)*[div(u*phi) + u.grad(phi)]
#     - theta3*div(gamma*eps*grad(phi)) + theta4*div(gamma*(phi - phi^2)*nhat) = 0
#
# phi is the phase field, 0 in one phase and 1 in the other, and nhat is the
# unit interface normal grad(phi)/|grad(phi)|. The diffusion and sharpening
# terms act against each other and hold the interface at a tanh profile of
# fixed width.
#
# Scheme: the same as Example 6 / Example 7. Forward Euler in time, so each
# step is one linear solve with the constant mass matrix; every spatial term is
# evaluated at the previous level, held in phio.
#
# Stabilization (derivation P29, "Stabilized Weak Form with Residual Based
# Closure"): three residual-based terms are added to the weak form,
#   + theta2/2*(u*tau*R, grad v) + theta2/2*(tau*R, div(u v))   along the flow
#   + theta4*(vs*tau*R, grad v)                                  along the sharpening velocity
# with vs = gamma*(1 - 2*phi)*nhat. Here div(u) = 0, so the three collapse to
#   + ((theta2*u + theta4*vs)*tau*R, grad v).
# Everything in them is taken at the previous level, as the rest of the explicit
# scheme is: the residual uses the lagged time derivative (phio - phioo)/dt, so
# the matrix stays the constant mass matrix. For linear elements the second
# derivatives and the curvature of the element-wise normal drop out, and the
# two linearized sharpening pieces of the derivation's L(phi) combine at
# phi = phis into gamma*(1 - 2*phio)*nhat.grad(phio):
#   R   = theta1*(phio - phioo)/dt + theta2*u.grad(phio)
#         + theta4*gamma*(1 - 2*phio)*nhat.grad(phio)
#   tau = ((theta1/dt)^2 + theta2^2 u.G.u + theta4^2 vs.G.vs
#          + C_I theta3^2 (gamma eps)^2 G:G)^(-1/2),  G = 4/h^2 I,  C_I = 36.
# nhat is a unit vector, so vs.G.vs = 4 gamma^2 (1 - 2 phio)^2 / h^2 and tau
# depends only on phio and position. It is therefore computed as a nodal field
# in Julia after every step and handed to Finch as a known variable, which
# avoids expanding it symbolically in the assembly (Example 13 showed that cost).
#
# Test problem: a disk of radius R with a rectangular slot cut into it, carried
# by the rotation u = (0.5 - y, x - 0.5) about the centre of the unit square.
# One revolution takes T = 2*pi, after which the exact solution is the initial
# field, so the error is measured against the initial condition. The disk never
# comes within 0.1 of a wall. The corners of the slot are not an equilibrium
# shape for a diffuse interface, so some rounding over a distance ~eps is
# expected; the body of the disk should return to where it started.
#
# Boundaries: the rotation crosses all four walls, but the disk stays 0.1 away
# from them and phi is ~1e-9 there. No condition is imposed on phi and no
# boundary integral is kept, which is the zero-total-flux (closed box) condition:
# the discrete mass is then conserved exactly. (Dirichlet phi = 0 instead would
# delete the tanh tails that touch the wall: measured -5e-3 of the mass on a
# 50 x 50 mesh and -1e-5 on 100 x 100, where the tails are longer.)
#
# Velocity: div(u) = 0 for solid-body rotation, so the div(u)*v part of the
# skew-symmetric form is dropped; the two halves of the advection term are
# otherwise written as in Example 6.
#
# Continuous Galerkin, linear elements on a quad mesh.
#==============================================================================#

using Finch

initFinch("cdi2dzalesakvms");

# results are written to output/, a folder kept in the repo
useLog("cdi2dzalesakvmslog", dir="output", level=3)

domain(2)
functionSpace(order=1)
timeStepper(EULER_EXPLICIT)

# The system matrix is the (symmetric positive definite) mass matrix, so
# plain conjugate gradients converge in a few iterations from the previous
# solution. The tight tolerance keeps the mass drift of 25,000 steps at round-off.
linAlgOptions(iterative=true, method="CG", pc="NONE", reltol=1e-12)

nel = 200
mesh(QUADMESH, elsperdim=nel)      # unit square

phi   = variable("phi")    # unknown at the new time level
phio  = variable("phio")   # known field: previous time level
phioo = variable("phioo")  # known field: two levels back (lagged time derivative)
tau   = variable("tau")    # known field: stabilization parameter, nodal
testSymbol("v")

# Physical parameters
dx   = 1.0/nel
umax = 0.5*sqrt(2)         # largest |u| in the square (at the corners)
gam  = 1.0                 # sharpening parameter gamma
eps  = dx                  # interface thickness (eps/dx >= (u/gamma + 1)/2)

# Slotted disk: centre, radius, slot width and the y-range of the slot. The
# slot runs from the bottom of the disk up past its centre.
xc  = 0.5
yc  = 0.75
R   = 0.15
sw  = 0.05                 # slot width
sy1 = yc - R               # slot bottom (the disk's bottom)
sy2 = yc + 0.10            # slot top

# Signed distance to the slotted disk, negative inside. The disk minus the slot
# is "inside disk AND outside slot", whose distance is the larger of the two
# signed distances (exact away from the corners, a good approximation at them).
function zalesak_sdf(x, y)
    d_disk = sqrt((x - xc)^2 + (y - yc)^2) - R           # negative inside the disk
    # signed distance to the slot rectangle [xc-sw/2, xc+sw/2] x [sy1, sy2]
    ddx = max(abs(x - xc) - sw/2, 0.0)
    ddy = max(max(sy1 - y, y - sy2), 0.0)
    outside = sqrt(ddx^2 + ddy^2)
    inside  = min(max(abs(x - xc) - sw/2, max(sy1 - y, y - sy2)), 0.0)
    d_slot  = outside + inside                            # negative inside the slot
    return max(d_disk, -d_slot)
end

# Initial condition: equilibrium tanh profile across the shape's boundary
zalesak_ic(x, y) = 0.5*(1 - tanh(zalesak_sdf(x, y)/(2*eps)))

# Zero total flux through the walls (NO_BC = nothing imposed, no boundary
# integral). Finch needs the declaration even though nothing is imposed.
boundary(phi, 1, NO_BC)

# Finch evaluates this string with x and y in scope
initial(phi, "Main.zalesak_ic(x, y)")

# Time stepping. Explicit diffusive limit in 2D (lumped mass):
# dt <= dx^2/(4*gamma*eps). Consistent mass matrix needs roughly 1/3 of that.
# One full revolution.
T  = 2*pi
dt = 0.2*dx*dx/(4*gam*eps)
nsteps = Int(round(T/dt))
dt = T/nsteps                      # land exactly on one revolution
setSteps(dt, nsteps)

# Coefficients (theta as in the derivation, dimensional form)
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", ["0.5 - y", "x - 0.5"], type=VECTOR)   # solid-body rotation
coefficient("gamma", gam)
coefficient("epsilon", eps)
coefficient("delta", 1e-12)                     # regularizes |grad(phi)| in the bulk
coefficient("idt", 1.0/dt)                      # 1/dt for the lagged time derivative

# The stabilization parameter at the nodes, from phio and the node position.
# Element metric on the uniform quad mesh with parent element [-1,1]^2:
# G = 4/h^2 I, u.G.u = 4|u|^2/h^2, vs.G.vs = 4 gamma^2 (1-2 phio)^2/h^2, G:G = 2 (4/h^2)^2.
CI = 36.0
Gm = 4.0/(dx*dx)
function update_tau!(tauvals, phiovals, nodes)
    for i = 1:length(tauvals)
        xi = nodes[1, i]; yi = nodes[2, i]
        u2 = (0.5 - yi)^2 + (xi - 0.5)^2
        vs2 = gam*gam*(1 - 2*phiovals[i])^2
        tauvals[i] = 1.0/sqrt((1.0/dt)^2 + u2*Gm + vs2*Gm + CI*(gam*eps)^2*2*Gm*Gm)
    end
end

# The unit interface normal, named here so the weak form below stays readable.
parameter("nhat", "grad(phio) ./ sqrt.(dot(grad(phio), grad(phio)) .+ delta)")

# Snapshots for the time series: 60 over the revolution plus t = 0.
snap_every = max(1, nsteps ÷ 60)
snapshots  = Vector{Vector{Float64}}()
snap_times = Float64[]
step_count = Ref(0)

# After each step the history shifts, tau is refreshed from the new previous
# level, and a snapshot is stored when due. Must be declared before weakForm.
@postStepFunction(
    begin
        phioo.values .= phio.values
        phio.values  .= phi.values
        update_tau!(tau.values, phio.values, Finch.finch_state.grid_data.allnodes)
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
#   - theta2/2*(phio, u.grad v)                       other half (div(u) = 0)
#   + theta3*(gamma*eps*grad phio, grad v)            diffusion
#   - theta4*(gamma*(phio - phio^2)*nhat, grad v)     sharpening
#   + ((theta2*u + theta4*vs)*tau*R, grad v)          VMS stabilization
# No boundary term and no boundary condition: zero total flux through the walls.
# The normal appears twice in the stabilization term (in vs and in R), and Finch
# rejects a parameter with a derivative used more than once, so it is written
# out inline there (NHAT below) instead of through the nhat parameter.
NHAT = "(grad(phio) ./ sqrt.(dot(grad(phio), grad(phio)) .+ delta))"
RES  = "(theta1*idt*(phio - phioo) + theta2*dot(u, grad(phio)) + theta4*gamma*(1 - 2*phio) .* dot($NHAT, grad(phio)))"
weakForm(phi, "theta1*Dt(phi)*v" *
              " - 0.5*theta2*phio*dot(u, grad(v))" *
              " - 0.5*theta2*phio*dot(u, grad(v))" *
              " + theta3*gamma*epsilon*dot(grad(phio), grad(v))" *
              " - theta4*gamma*(phio - phio*phio) .* dot(nhat, grad(v))" *
              " + theta2*tau*$RES .* dot(u, grad(v))" *
              " + theta4*gamma*(1 - 2*phio)*tau*$RES .* dot($NHAT, grad(v))")

# Evaluate the initial condition now so it can be stored as the t = 0 snapshot
evalInitialConditions()
phio.values  .= phi.values
phioo.values .= phi.values         # zero lagged time derivative on the first step
phi0 = copy(phi.values[:])         # the exact solution after one revolution
update_tau!(tau.values, phio.values, Finch.finch_state.grid_data.allnodes)
push!(snapshots, copy(phi.values[:]))
push!(snap_times, 0.0)

# Lumped finite element weights (half at an edge node, quarter at a corner)
nodes = Finch.finch_state.grid_data.allnodes
wts = fill(dx*dx, size(nodes, 2))
for i = 1:size(nodes, 2)
    if nodes[1,i] < 1e-12 || nodes[1,i] > 1 - 1e-12; wts[i] *= 0.5; end
    if nodes[2,i] < 1e-12 || nodes[2,i] > 1 - 1e-12; wts[i] *= 0.5; end
end
mass0 = sum(wts .* phi0)

solve(phi)

finalizeFinch()

# Error after one revolution, against the initial field
x = nodes[1, :]
y = nodes[2, :]
phivals = phi.values[:]
err = phivals .- phi0

maxerr = maximum(abs.(err))
l1err  = sum(wts .* abs.(err))            # the usual Zalesak measure
massT  = sum(wts .* phivals)
println("dt = " * string(dt) * ", steps = " * string(nsteps) * " (one revolution, T = 2*pi)")
println("max error = " * string(maxerr))
println("L1 error = " * string(l1err))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))
println("relative mass change = " * string((massT - mass0)/mass0))

using DelimitedFiles

data = [x y phivals phi0]
writedlm("output/solution.csv", data, ',')

# Time series: columns 1-2 = x, y; columns 3.. = phi at each snapshot time
ts = hcat(x, y, snapshots...)
writedlm("output/solution_timeseries.csv", ts, ',')
writedlm("output/snapshot_times.csv", snap_times, ',')
println(string(length(snap_times)) * " snapshots written")
