#==============================================================================#
# 2D Conservative Diffusive Interface (CDI) model, semi-implicit: BDF2 in time
# (BDF1 on the first step) with the sharpening term linearized and two Picard
# iterations per step, each one linear solve. Finch V2 (branch dev); the scheme
# is in scheme.jl (shared with sharpening_sweep.jl), the shared helpers and the
# FINCH_TARGET switch (CPU Julia, or MLIR on the CPU/GPU) in ../common.jl.
#
#   theta1*phi_t + theta2*div(u*phi) - theta3*div(gamma*eps*grad(phi))
#     + theta4*div(gamma*(phi - phi^2)*nhat) = 0,   nhat = grad(phi)/|grad(phi)|
#
# The normal is recovered at the nodes by an L2 projection of the gradient of
# the lagged field (consistent mass matrix). Run with ./run.sh; CDI_CASE picks
# the benchmark (see scheme.jl):
#   circle     (default) a circle carried through the square by u = (1, 0),
#              compared with the translated profile at T = 0.4 (mesh 300 x 300)
#   zalesak    Zalesak's slotted disk, one revolution, compared with the initial
#              field: max and L1 (lumped) error (mesh 200 x 200)
#   sheardrop  a drop in the time-reversed vortex of Mirjalili et al. 2020,
#              one period, the same errors plus the most stretched state at
#              T/2 (mesh 200 x 200; paper512.sh runs the paper's parameters)
# CDI_NEL overrides the mesh, CDI_OUTDIR the output folder (default
# output/<case>/), CDI_GAMMA and CDI_EPS_DX the interface parameters.
#
# Output: solution.csv (x, y, phi at T, the reference field, and phi at T/2
# for the shear drop), solution_timeseries.csv, snapshot_times.csv.
#==============================================================================#

include("scheme.jl")

const CASE = get(ENV, "CDI_CASE", "circle")
CASE in CASES || error("CDI_CASE must be circle, zalesak or sheardrop, got \"$CASE\"")
const CFL = 0.25
const OUT_DIR = joinpath(@__DIR__, get(ENV, "CDI_OUTDIR", joinpath("output", CASE)))
mkpath(OUT_DIR)

r = run_case(CFL, 2; case = CASE, nel = parse(Int, get(ENV, "CDI_NEL", string(case_nel(CASE)))),
             proj_dir = joinpath(OUT_DIR, "finch"), verbose = true)
r.ok || error("the run failed")
x, y, phivals, T = r.x, r.y, r.phi, r.T

if CASE == "circle"
    # Error against the translated exact profile
    xc = XC0 + 1.0*T
    yc = YC0 + 0.0*T
    rr = sqrt.((x .- xc).^2 .+ (y .- yc).^2)
    exact = 0.5 .* (1 .- tanh.((rr .- R) ./ (2*r.eps)))
    maxerr = maximum(abs.(phivals .- exact))
    println("dt = " * string(r.dt) * ", steps = " * string(r.nsteps) * " (CFL " * string(CFL) * ")")
    println("max error = " * string(maxerr))
    println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))
    println("relative mass change = " * string((r.mass[end] - r.mass[1])/r.mass[1]))
    data = [x y phivals exact]
else
    # Error after one revolution or one period, against the initial field
    phi0 = r.snapshots[1]
    dx   = 1.0/round(Int, 1/minimum(diff(sort(unique(x)))))
    wts  = lumped_weights(x, y, dx)      # half at an edge node, quarter at a corner
    err  = phivals .- phi0
    maxerr = maximum(abs.(err))
    l1err  = sum(wts .* abs.(err))      # the usual Zalesak measure / the shape error of the paper
    if CASE == "zalesak"
        println("dt = " * string(r.dt) * ", steps = " * string(r.nsteps) * " (CFL " * string(CFL) * ", one revolution)")
        println("max error = " * string(maxerr))
        println("L1 error = " * string(l1err))
        data = [x y phivals phi0]
    else
        nel = round(Int, 1/dx)
        phihalf = r.snapshots[argmin(abs.(r.snap_times .- T/2))]   # the most stretched state
        println("mesh " * string(nel) * " x " * string(nel) * ", eps/dx = " * string(r.eps/dx) * ", gamma/umax = " * string(r.gam/r.umax))
        println("dt = " * string(r.dt) * ", steps = " * string(r.nsteps) * " (CFL " * string(CFL) * ", one period T = " * string(T) * ")")
        println("max error = " * string(maxerr))
        println("L1 (shape) error = " * string(l1err))
        data = [x y phivals phi0 phihalf]
    end
    println("min phi = " * string(minimum(phivals)) * ", max phi = " * string(maximum(phivals)))
    println("relative mass change = " * string((r.mass[end] - r.mass[1])/r.mass[1]))
end

writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')
# Time series: columns 1-2 = x, y; columns 3.. = phi at each snapshot time
write_timeseries(OUT_DIR, "solution", [x y], reduce(hcat, r.snapshots))
write_snapshot_times(OUT_DIR, r.snap_times)
println(string(length(r.snap_times)) * " snapshots written")
