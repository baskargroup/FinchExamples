#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, solved explicitly in Finch:
# Zalesak's slotted disk in solid-body rotation.
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
# and refreshed after every step (Examples 5 and 6 do the same). Reference Examples 10-16
# take the normal element by element from the gradient of the linear basis;
# that normal is constant inside each element and jumps between elements, and
# on this case it produced grid-scale wiggles along the slot walls, shed
# fragments and a widened slot (L1 error 1.2e-2, slot 0.059 wide instead of
# 0.050; elementwise_normal.sh reproduces it). With the nodal normal the
# walls are smooth, nothing is shed and the slot keeps its width (L1 error
# 3.4e-3). Mirjalili et al. (JCP 2020) likewise use a continuous normal, from
# central differences at cell centres.
#
# Scheme: the same as Example 6. Forward Euler in time, so each step is one
# linear solve with the constant mass matrix; every spatial term is evaluated
# at the previous level, held in phio.
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

initFinch("cdi2dzalesak");

# Results are written to output/, a folder kept in the repo. The mesh, the
# interface parameters and the output folder can be overridden from the
# environment (see paper512.sh); the defaults are the 200 x 200 case.
outdir = get(ENV, "CDI_OUTDIR", "output")
mkpath(outdir)
useLog("cdi2dzalesaklog", dir=outdir, level=3)

domain(2)
functionSpace(order=1)
timeStepper(EULER_EXPLICIT)

# The system matrix is the (symmetric positive definite) mass matrix, so
# plain conjugate gradients converge in a few iterations from the previous
# solution. The tight tolerance keeps the mass drift of 25,000 steps at round-off.
linAlgOptions(iterative=true, method="CG", pc="NONE", reltol=1e-12)

nel = parse(Int, get(ENV, "CDI_NEL", "200"))
mesh(QUADMESH, elsperdim=nel)      # unit square

phi  = variable("phi")     # unknown at the new time level
phio = variable("phio")    # known field: previous time level
nhat = variable("nhat", type=VECTOR)   # known field: nodal unit normal, built from phio
testSymbol("v")

# Physical parameters
dx   = 1.0/nel
umax = 0.5*sqrt(2)         # largest |u| in the square (at the corners)
# sharpening parameter gamma: 1 by default, CDI_GAMMA * umax when given
gam  = haskey(ENV, "CDI_GAMMA") ? parse(Float64, ENV["CDI_GAMMA"])*umax : 1.0
eps  = parse(Float64, get(ENV, "CDI_EPS_DX", "1.0"))*dx   # interface thickness (eps/dx >= (u/gamma + 1)/2)

# Slotted disk: centre, radius, slot width and the y-range of the slot. The
# slot opens at the bottom of the disk (y = 0.6) and runs up past its centre.
# The slot RECTANGLE used in the distance function is extended 0.05 below the
# disk: it removes nothing more from the disk, but if its bottom edge sat at
# the disk's bottom the distance function would see a boundary there and the
# tanh profile would close the slot mouth with a ridge of phi = 0.5.
xc  = 0.5
yc  = 0.75
R   = 0.15
sw  = 0.05                 # slot width
sy1 = yc - R - 0.05        # slot rectangle bottom, below the disk (see above)
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
# (the 1e-12 that regularizes |grad(phi)| in the bulk is passed to update_nhat! below)

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

# Snapshots for the time series: 60 over the revolution plus t = 0.
snap_every = max(1, nsteps ÷ 60)
snapshots  = Vector{Vector{Float64}}()
snap_times = Float64[]
step_count = Ref(0)

# After each step the new solution becomes the "previous" field, the nodal
# normal is rebuilt from it, and a snapshot is stored when due. Must be
# declared before weakForm.
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
#   - theta2/2*(phio, u.grad v)                       other half (div(u) = 0)
#   + theta3*(gamma*eps*grad phio, grad v)            diffusion
#   - theta4*(gamma*(phio - phio^2)*nhat, grad v)     sharpening, nodal nhat
# No boundary term and no boundary condition: zero total flux through the walls.
weakForm(phi, "theta1*Dt(phi)*v" *
              " - 0.5*theta2*phio*dot(u, grad(v))" *
              " - 0.5*theta2*phio*dot(u, grad(v))" *
              " + theta3*gamma*epsilon*dot(grad(phio), grad(v))" *
              " - theta4*gamma*(phio - phio*phio) .* dot(nhat, grad(v))")

# Evaluate the initial condition now so it can be stored as the t = 0 snapshot
evalInitialConditions()
phio.values .= phi.values
recover_nhat!(nhat.values, phio.values, 1e-12)
phi0 = copy(phi.values[:])         # the exact solution after one revolution
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
writedlm(joinpath(outdir, "solution.csv"), data, ',')

# Time series: columns 1-2 = x, y; columns 3.. = phi at each snapshot time
ts = hcat(x, y, snapshots...)
writedlm(joinpath(outdir, "solution_timeseries.csv"), ts, ',')
writedlm(joinpath(outdir, "snapshot_times.csv"), snap_times, ',')
println(string(length(snap_times)) * " snapshots written")
