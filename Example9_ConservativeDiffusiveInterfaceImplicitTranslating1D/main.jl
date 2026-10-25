#==============================================================================#
# 1D Conservative Diffusive Interface (CDI) model, semi-implicit, with
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
# element gradient of the Picard iterate phis, (g, v) = (phis_x, v), refreshed
# before each of the two solves; nhat = g/sqrt(g^2 + delta). Finch V2 has no
# lumped mass matrix, so the projection uses the consistent one (the V1
# version of this example lumped it). The element-wise normal from the
# gradient of the linear basis is constant per element and jumps between
# elements; in 2D it grew a grid-scale zigzag on grid-aligned fronts
# (Example 10), and the nodal normal removes it. In 1D the two agree.
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
# as a step_schedule: phi is solved as a quasi-static equation (no Dt) whose
# time derivative is the BDF formula over the known fields phio and phioo,
# and the history, the lagged field and the BDF coefficients are moved along
# by update entries:
#   phis <- phio, project g, solve phi       iteration 1
#   phis <- phi,  project g, solve phi       iteration 2
#   phioo <- phio, phio <- phi, BDF1 -> BDF2
#
# Test problem: one interface carried by a uniform velocity. The equilibrium
# profile translates without changing shape, so the exact solution is
#   phi(x,t) = 0.5*(1 + tanh((x - x0 - u*t)/(2*eps)))
#
# Continuous Galerkin, linear elements.
#==============================================================================#

using DelimitedFiles

include("scheme.jl")

# Time step. Diffusion is implicit, so the explicit limit dx^2/(2*gamma*eps)
# does not apply; dt is set by accuracy through the advective CFL number.
# CFL 0.25 keeps phi within [0,1] to round-off and sits near the spatial error
# floor; CFL 0.5 still beats the explicit scheme but overshoots by 4e-4, and
# CFL 1 and above lose accuracy quickly (see notes/).
const CFL = 0.25

const OUT_DIR = joinpath(@__DIR__, "output")
r = run_case(CFL, 2; proj_dir = joinpath(OUT_DIR, "finch"), snapshots = 20)

# Compare with the exact profile at the final time
x, phivals = r.x, r.phi
exact = 0.5 .* (1 .+ tanh.((x .- X0 .- UVEL*T) ./ (2*EPS)))
maxerr = maximum(abs.(phivals .- exact))
println("dt = " * string(r.dt) * ", steps = " * string(r.nsteps) * " (CFL " * string(CFL) * ")")
println("max error = " * string(maxerr))
println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))

data = [x phivals exact]
writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')

# Time series: column 1 is x, the rest is phi at each snapshot time
writedlm(joinpath(OUT_DIR, "solution_timeseries.csv"), hcat(x, r.snapshots...), ',')
writedlm(joinpath(OUT_DIR, "snapshot_times.csv"), r.snap_times, ',')
println(string(length(r.snap_times)) * " snapshots written")
