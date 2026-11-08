#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, semi-implicit, with
# Finch V2 (Julia target).
#
#   theta1*phi_t + theta2*(1/2)*[div(u*phi) + u.grad(phi)]
#     - theta3*div(gamma*eps*grad(phi)) + theta4*div(gamma*(phi - phi^2)*nhat) = 0
#
# phi is the phase field, 0 in one phase and 1 in the other, and nhat is the
# unit interface normal grad(phi)/|grad(phi)|. The diffusion and sharpening
# terms act against each other and hold the interface at a tanh profile of
# fixed width.
#
# Interface normal: recovered at the NODES by an L2 projection of the finite
# element gradient of the Picard iterate phis, one field per component,
# refreshed before each of the two solves (consistent mass; Finch V2 has no
# lumped one, the V1 version of this example lumped it). With the element-wise
# normal from the gradient of the linear basis (constant per element, jumping
# between elements) this case grew a 3-element zigzag along the leading edge
# from t ~ 0.3 on this mesh, max error 0.128, with none of the time-scheme or
# linearization changes curing it; with the nodal normal the front is flat to
# 0.02 elements and the max error is 0.0068.
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
# Finch V2 has forward and backward Euler only, so one time step is written
# as a step_schedule, as in Example 9 (the scheme is in scheme.jl, shared with
# sharpening_sweep.jl). The matrix is not symmetric because of advection, so
# BiCGSTAB is used.
#
# Test problem: a circle of radius R carried by a uniform velocity, the same
# case as Example 6. The profile translates with almost no change of shape, so
#   phi(x,y,t) = 0.5*(1 - tanh((r - R)/(2*eps)))
# measured from the moving centre is close to the solution; a curved interface
# is not an equilibrium profile, so it is not exact.
#
# Boundaries: phi = 0 is imposed at the inflow x = 0. At the outflow x = 1 the
# convective boundary flux theta2*(u.n)*phi is kept and the diffusive and
# sharpening fluxes are dropped (the do-nothing outflow condition), so mass
# leaves at the rate the flow carries it out. The V1 version evaluated that
# flux from the lagged field phis; here it is implicit (scheme.jl). Top and
# bottom have u.n = 0.
#
# Continuous Galerkin, linear elements on a quad mesh.
#==============================================================================#

using DelimitedFiles

include("scheme.jl")

# Time step. Diffusion is implicit, so the explicit limit dx^2/(4*gamma*eps)
# does not apply; dt is set by accuracy through the advective CFL number.
# CFL 0.25 kept phi within [0,1] and sat near the spatial error floor in the
# 1D study (Example 9).
const CFL = 0.25

const OUT_DIR = joinpath(@__DIR__, "output")
r = run_case(CFL, 2; nel = parse(Int, get(ENV, "CDI_NEL", string(NEL))),
             proj_dir = joinpath(OUT_DIR, "finch"), snapshots = 60, verbose = true)

# Error check against the translated exact profile
x, y, phivals = r.x, r.y, r.phi
eps = 1.0/round(Int, 1/minimum(diff(sort(unique(x)))))
xc = XC0 + UX*T
yc = YC0 + UY*T
rr = sqrt.((x .- xc).^2 .+ (y .- yc).^2)
exact = 0.5 .* (1 .- tanh.((rr .- R) ./ (2*eps)))

maxerr = maximum(abs.(phivals .- exact))
println("dt = " * string(r.dt) * ", steps = " * string(r.nsteps) * " (CFL " * string(CFL) * ")")
println("max error = " * string(maxerr))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))
println("relative mass change = " * string((r.mass[end] - r.mass[1])/r.mass[1]))

data = [x y phivals exact]
writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')

# Time series: columns 1-2 = x, y; columns 3.. = phi at each snapshot time
writedlm(joinpath(OUT_DIR, "solution_timeseries.csv"), hcat(x, y, r.snapshots...), ',')
writedlm(joinpath(OUT_DIR, "snapshot_times.csv"), r.snap_times, ',')
println(string(length(r.snap_times)) * " snapshots written")
