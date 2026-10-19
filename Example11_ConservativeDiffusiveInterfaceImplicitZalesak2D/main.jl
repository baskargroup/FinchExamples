#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, semi-implicit, in Finch:
# Zalesak's slotted disk in solid-body rotation.
#
#   theta1*phi_t + theta2*(1/2)*[div(u*phi) + u.grad(phi)]
#     - theta3*div(gamma*eps*grad(phi)) + theta4*div(gamma*(phi - phi^2)*nhat) = 0
#
# The scheme of Example 10 on the problem of Example 7.
#
# Time discretization (semi-implicit, two Picard iterations per step):
#   - BDF2 for the time derivative, BDF1 on the first step:
#       phi_t ~ (a0*phi + a1*phio + a2*phioo)/dt
#   - advection and diffusion implicit (linear in phi)
#   - sharpening linearized with a lagged field phis, as derived:
#       gamma*(phi - phi^2)*nhat ~ gamma*phi*(1 - phis)*nhat(phis)
#   - iteration 1 uses phis = phio, iteration 2 uses phis = the result of
#     iteration 1. Each iteration is one linear solve (GMRES); the matrix
#     depends on phis, so it is rebuilt for iteration 2.
#
# Interface normal: recovered at the NODES by a lumped L2 projection of the
# finite element gradient of phis, refreshed before each solve (as in
# Examples 5-10).
#
# Test problem (as Example 7): a disk of radius R with a rectangular slot,
# carried by the rotation u = (0.5 - y, x - 0.5) about the centre of the unit
# square. After one revolution, T = 2*pi, the exact solution is the initial
# field, so the error is measured against it.
#
# Boundaries (as Example 7): no condition on phi and no boundary integral,
# the zero-total-flux (closed box) condition; the disk stays 0.1 away from
# the walls and the discrete mass is conserved.
#
# Velocity: div(u) = 0 for solid-body rotation, so the div(u)*v part of the
# skew-symmetric advection form is dropped.
#
# Continuous Galerkin, linear elements on a quad mesh.
#==============================================================================#

using Finch

initFinch("cdi2dzalesakimplicit");

# results are written to output/, a folder kept in the repo
useLog("cdi2dzalesakimplicitlog", dir="output", level=2)

domain(2)
functionSpace(order=1)
# No timeStepper() call: Finch sees a steady problem. Time stepping is below.

# The matrix changes every Picard iteration and is not symmetric (advection),
# so it is assembled each time and solved with GMRES (see Example 10).
linAlgOptions(iterative=true, method="GMRES", pc="NONE", reltol=1e-12)

nel = 200
mesh(QUADMESH, elsperdim=nel)      # unit square

phi   = variable("phi")    # unknown: the current Picard iterate
phio  = variable("phio")   # known: previous time level
phioo = variable("phioo")  # known: two time levels back
phis  = variable("phis")   # known: lagged field for the sharpening term
nhat  = variable("nhat", type=VECTOR)   # known: nodal unit normal, built from phis
testSymbol("v")

# Physical parameters (as Example 7)
dx   = 1.0/nel
umax = 0.5*sqrt(2)         # largest |u| in the square (at the corners)
gam  = 1.0                 # sharpening parameter gamma
eps  = dx                  # interface thickness (eps/dx >= (u/gamma + 1)/2)

# Slotted disk (as Example 7; the slot rectangle extends 0.05 below the disk
# so that the slot mouth is open)
xc  = 0.5
yc  = 0.75
R   = 0.15
sw  = 0.05                 # slot width
sy1 = yc - R - 0.05        # slot rectangle bottom, below the disk
sy2 = yc + 0.10            # slot top

# Signed distance to the slotted disk, negative inside
function zalesak_sdf(x, y)
    d_disk = sqrt((x - xc)^2 + (y - yc)^2) - R
    ddx = max(abs(x - xc) - sw/2, 0.0)
    ddy = max(max(sy1 - y, y - sy2), 0.0)
    outside = sqrt(ddx^2 + ddy^2)
    inside  = min(max(abs(x - xc) - sw/2, max(sy1 - y, y - sy2)), 0.0)
    d_slot  = outside + inside
    return max(d_disk, -d_slot)
end

zalesak_ic(x, y) = 0.5*(1 - tanh(zalesak_sdf(x, y)/(2*eps)))

# Zero total flux through the walls. Finch needs the declaration even though
# nothing is imposed.
boundary(phi, 1, NO_BC)

initial(phi, "Main.zalesak_ic(x, y)")

# Time step from the advective CFL number with the largest speed in the
# square; CFL 0.25 as in Examples 9 and 10. One full revolution.
T   = 2*pi
CFL = 0.25
dt  = CFL*dx/umax
nsteps = Int(round(T/dt))
dt  = T/nsteps                     # land exactly on one revolution

# BDF coefficients, divided by dt, read by Finch at every assembly
a0_dt = 1.0/dt;  a1_dt = -1.0/dt;  a2_dt = 0.0           # BDF1
bdf2!() = (global a0_dt = 1.5/dt; global a1_dt = -2.0/dt; global a2_dt = 0.5/dt)

# Coefficients (theta as in the derivation, dimensional form)
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", ["0.5 - y", "x - 0.5"], type=VECTOR)   # solid-body rotation
coefficient("gamma", gam)
coefficient("epsilon", eps)
coefficient("a0dt", "Main.a0_dt")
coefficient("a1dt", "Main.a1_dt")
coefficient("a2dt", "Main.a2_dt")

# Nodal interface normal by finite element gradient recovery (lumped L2
# projection), identical to Examples 7 and 10
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
        J = g_geo.J[e]                      # constant Jacobian on this mesh
        rx = J.rx[1]; ry = J.ry[1]; sx = J.sx[1]; sy = J.sy[1]
        detj = g_geo.detJ[1, e]
        for i = 1:g_Np
            g_pe[i] = p[g_grid.loc2glb[i, e]]
        end
        for q = 1:g_Nqp
            dr = 0.0; ds = 0.0
            for i = 1:g_Np
                dr += g_refel.Qr[q, i] * g_pe[i]
                ds += g_refel.Qs[q, i] * g_pe[i]
            end
            g_dr[q] = dr; g_ds[q] = ds
        end
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

# Weak form of one Picard iteration (residual = 0):
#   theta1*(a0/dt*phi, v) + theta1*((a1*phio + a2*phioo)/dt, v)   time derivative
#   - theta2/2*(u*phi, grad v)                                    advection, conservative half
#   - theta2/2*(phi, u.grad v)                                    other half (div(u) = 0)
#   + theta3*(gamma*eps*grad phi, grad v)                         diffusion
#   - theta4*(gamma*phi*(1 - phis)*nhat(phis), grad v)            sharpening, lagged
# No boundary term and no boundary condition: zero total flux through the walls.
weakForm(phi, "theta1*a0dt*phi*v" *
              " + theta1*(a1dt*phio + a2dt*phioo)*v" *
              " - 0.5*theta2*phi*dot(u, grad(v))" *
              " - 0.5*theta2*phi*dot(u, grad(v))" *
              " + theta3*gamma*epsilon*dot(grad(phi), grad(v))" *
              " - theta4*gamma*phi*(1 - phis) .* dot(nhat, grad(v))")

# Snapshots for the time series: 60 over the revolution plus t = 0
snap_every = max(1, nsteps ÷ 60)
snapshots  = Vector{Vector{Float64}}()
snap_times = Float64[]

# Start: all history levels equal the initial condition
evalInitialConditions()
phio.values  .= phi.values
phioo.values .= phi.values
phis.values  .= phi.values
phi0 = copy(phi.values[:])         # the exact solution after one revolution
push!(snapshots, copy(phi0))
push!(snap_times, 0.0)

# Lumped finite element weights (half at an edge node, quarter at a corner)
nodes = Finch.finch_state.grid_data.allnodes
wts = fill(dx*dx, size(nodes, 2))
for i = 1:size(nodes, 2)
    if nodes[1,i] < 1e-12 || nodes[1,i] > 1 - 1e-12; wts[i] *= 0.5; end
    if nodes[2,i] < 1e-12 || nodes[2,i] > 1 - 1e-12; wts[i] *= 0.5; end
end
mass0 = sum(wts .* phi0)

# Time loop. Two Picard iterations per step, then shift the history.
t0 = time()
for step = 1:nsteps
    if step == 2
        bdf2!()                          # BDF1 only on the first step
    end
    phis.values .= phio.values           # iteration 1: lag with the previous level
    recover_nhat!(nhat.values, phis.values, 1e-12)
    solve(phi)
    phis.values .= phi.values            # iteration 2: lag with iteration 1
    recover_nhat!(nhat.values, phis.values, 1e-12)
    solve(phi)
    phioo.values .= phio.values          # shift history
    phio.values  .= phi.values
    if step % snap_every == 0
        push!(snapshots, copy(phi.values[:]))
        push!(snap_times, step*dt)
    end
    if step % 100 == 0 || step == nsteps
        println("step $step / $nsteps   t = $(round(step*dt, digits=4))   elapsed $(round(time() - t0, digits=1)) s")
        flush(stdout)
    end
end

finalizeFinch()

# Error after one revolution, against the initial field
x = nodes[1, :]
y = nodes[2, :]
phivals = phi.values[:]
err = phivals .- phi0

maxerr = maximum(abs.(err))
l1err  = sum(wts .* abs.(err))            # the usual Zalesak measure
massT  = sum(wts .* phivals)
println("dt = " * string(dt) * ", steps = " * string(nsteps) * " (CFL " * string(CFL) * ", one revolution)")
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
