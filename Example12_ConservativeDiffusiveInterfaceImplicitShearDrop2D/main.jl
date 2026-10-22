#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, semi-implicit, in Finch:
# a drop in a time-reversed shear (vortex) flow.
#
#   theta1*phi_t + theta2*(1/2)*[div(u*phi) + u.grad(phi)]
#     - theta3*div(gamma*eps*grad(phi)) + theta4*div(gamma*(phi - phi^2)*nhat) = 0
#
# The scheme of Examples 10 and 11 on the drop in shear flow test of
# Mirjalili et al., JCP 401 (2020) 109006, Section 4.3.1.
#
# Time discretization (semi-implicit, two Picard iterations per step):
#   - BDF2 for the time derivative, BDF1 on the first step:
#       phi_t ~ (a0*phi + a1*phio + a2*phioo)/dt
#   - advection and diffusion implicit (linear in phi); the velocity is
#     evaluated at the new time level t^{n+1}
#   - sharpening linearized with a lagged field phis, as derived:
#       gamma*(phi - phi^2)*nhat ~ gamma*phi*(1 - phis)*nhat(phis)
#   - iteration 1 uses phis = phio, iteration 2 uses phis = the result of
#     iteration 1. Each iteration is one linear solve (GMRES); the matrix
#     depends on phis, so it is rebuilt for iteration 2.
#
# Interface normal: recovered at the NODES by a lumped L2 projection of the
# finite element gradient of phis, refreshed before each solve (as in
# Examples 5-11).
#
# Test problem (Mirjalili et al. 2020, Section 4.3.1): a drop of radius
# R = 0.15 centred at (0.5, 0.75) in the unit square, in the velocity of the
# stream function
#     Psi(x, y, t) = (1/pi) sin^2(pi x) sin^2(pi y) cos(pi t/T),   T = 4,
#     u = -dPsi/dy = -sin^2(pi x) sin(2 pi y) cos(pi t/T)
#     v =  dPsi/dx =  sin(2 pi x) sin^2(pi y) cos(pi t/T)
# The vortex stretches the drop into a thin spiral up to t = T/2, then the
# flow reverses and undoes the stretching. At t = T the exact solution is the
# initial field, so the error is measured against it (the shape error of the
# paper, Eq. 41, is the L1 error below).
#
# Boundaries: no condition on phi and no boundary integral, the zero-total-
# flux (closed box) condition, as in Examples 7 and 11. The normal velocity
# vanishes on the walls, and the discrete mass is conserved.
#
# Velocity: div(u) = 0 (it comes from a stream function), so the div(u)*v
# part of the skew-symmetric advection form is dropped.
#
# Continuous Galerkin, linear elements on a quad mesh.
#==============================================================================#

using Finch

initFinch("cdi2dsheardropimplicit");

# Results are written to output/, a folder kept in the repo. The mesh, the
# interface parameters and the output folder can be overridden from the
# environment (see paper512.sh); the defaults are the 200 x 200 case.
outdir = get(ENV, "CDI_OUTDIR", "output")
mkpath(outdir)
useLog("cdi2dsheardropimplicitlog", dir=outdir, level=2)

domain(2)
functionSpace(order=1)
# No timeStepper() call: Finch sees a steady problem. Time stepping is below.

# The matrix changes every Picard iteration and is not symmetric (advection),
# so it is assembled each time and solved with GMRES (see Example 10).
linAlgOptions(iterative=true, method="GMRES", pc="NONE", reltol=1e-12)

nel = parse(Int, get(ENV, "CDI_NEL", "200"))
mesh(QUADMESH, elsperdim=nel)      # unit square

phi   = variable("phi")    # unknown: the current Picard iterate
phio  = variable("phio")   # known: previous time level
phioo = variable("phioo")  # known: two time levels back
phis  = variable("phis")   # known: lagged field for the sharpening term
nhat  = variable("nhat", type=VECTOR)   # known: nodal unit normal, built from phis
testSymbol("v")

# Physical parameters
dx   = 1.0/nel
umax = 1.0                 # largest |u| over space and time (at t = 0 and T)
gam  = parse(Float64, get(ENV, "CDI_GAMMA", "1.0"))*umax    # sharpening parameter gamma
eps  = parse(Float64, get(ENV, "CDI_EPS_DX", "1.0"))*dx   # interface thickness (eps/dx >= (u/gamma + 1)/2)

# Drop (Mirjalili et al. 2020, Section 4.3.1)
xc = 0.5
yc = 0.75
R  = 0.15

drop_ic(x, y) = 0.5*(1 - tanh((sqrt((x - xc)^2 + (y - yc)^2) - R)/(2*eps)))

# Zero total flux through the walls. Finch needs the declaration even though
# nothing is imposed.
boundary(phi, 1, NO_BC)

initial(phi, "Main.drop_ic(x, y)")

# Time step from the advective CFL number with the largest speed; CFL 0.25
# as in Examples 9-11. One full period of the flow.
T   = 4.0
CFL = 0.25
dt  = CFL*dx/umax
nsteps = Int(round(T/dt))
dt  = T/nsteps                     # land exactly on t = T

# Time factor of the velocity, cos(pi*t/T), set to its value at the new time
# level before each step and read by Finch at every assembly
const VEL_C = Ref(1.0)
set_velocity!(t) = (VEL_C[] = cos(pi*t/T))
# The velocity is given to Finch as calls to these functions. Reading the time
# factor from a typed const Ref keeps them type-stable; an untyped global read
# from the coefficient strings made each step 1.3 times slower here (2.3
# times in the explicit Example 8).
ux(x, y) = -VEL_C[]*sin(pi*x)^2*sin(2*pi*y)
uy(x, y) =  VEL_C[]*sin(2*pi*x)*sin(pi*y)^2

# BDF coefficients, divided by dt, read by Finch at every assembly
a0_dt = 1.0/dt;  a1_dt = -1.0/dt;  a2_dt = 0.0           # BDF1
bdf2!() = (global a0_dt = 1.5/dt; global a1_dt = -2.0/dt; global a2_dt = 0.5/dt)

# Coefficients (theta as in the derivation, dimensional form)
coefficient("theta1", 1.0)
coefficient("theta2", 1.0)
coefficient("theta3", 1.0)
coefficient("theta4", 1.0)
coefficient("u", ["Main.ux(x, y)", "Main.uy(x, y)"], type=VECTOR)   # vortex
coefficient("gamma", gam)
coefficient("epsilon", eps)
coefficient("a0dt", "Main.a0_dt")
coefficient("a1dt", "Main.a1_dt")
coefficient("a2dt", "Main.a2_dt")

# Nodal interface normal by finite element gradient recovery (lumped L2
# projection), identical to Examples 7, 10 and 11
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

# Snapshots for the time series: 64 over the period plus t = 0
snap_every = max(1, nsteps ÷ 64)
snapshots  = Vector{Vector{Float64}}()
snap_times = Float64[]

# Start: all history levels equal the initial condition
evalInitialConditions()
phio.values  .= phi.values
phioo.values .= phi.values
phis.values  .= phi.values
phi0 = copy(phi.values[:])         # the exact solution at t = T
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
    set_velocity!(step*dt)               # velocity at the new time level
    phis.values .= phio.values           # iteration 1: lag with the previous level
    recover_nhat!(nhat.values, phis.values, 1e-12)
    solve(phi)
    phis.values .= phi.values            # iteration 2: lag with iteration 1
    recover_nhat!(nhat.values, phis.values, 1e-12)
    solve(phi)
    phioo.values .= phio.values          # shift history
    phio.values  .= phi.values
    if step == nsteps ÷ 2
        global phihalf = copy(phi.values[:])   # the most stretched state, t = T/2
    end
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

# Error after one period, against the initial field
x = nodes[1, :]
y = nodes[2, :]
phivals = phi.values[:]
err = phivals .- phi0

maxerr = maximum(abs.(err))
l1err  = sum(wts .* abs.(err))            # shape error, Eq. 41 of the paper
massT  = sum(wts .* phivals)
println("dt = " * string(dt) * ", steps = " * string(nsteps) * " (CFL " * string(CFL) * ", one period T = " * string(T) * ")")
println("max error = " * string(maxerr))
println("L1 (shape) error = " * string(l1err))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))
println("relative mass change = " * string((massT - mass0)/mass0))

using DelimitedFiles

data = [x y phivals phi0 phihalf]
writedlm(joinpath(outdir, "solution.csv"), data, ',')

# Time series: columns 1-2 = x, y; columns 3.. = phi at each snapshot time
ts = hcat(x, y, snapshots...)
writedlm(joinpath(outdir, "solution_timeseries.csv"), ts, ',')
writedlm(joinpath(outdir, "snapshot_times.csv"), snap_times, ',')
println(string(length(snap_times)) * " snapshots written")
