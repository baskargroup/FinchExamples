#==============================================================================#
# Temporal convergence study for the semi-implicit 1D CDI scheme of main.jl
# (Finch V2, the scheme in scheme.jl).
#
# The error against the exact tanh contains a spatial part (~5e-3 on this mesh)
# that does not shrink with dt, so it cannot show the order in time. Here each
# run is instead compared with a reference solution on the SAME mesh at a very
# small dt: the spatial error is identical in both and cancels, and what is
# left is the time-stepping error alone.
#
# Two variants are measured: the scheme as derived (two Picard iterations per
# step, expected second order) and a single-iteration variant (lagged field =
# previous level, expected first order).
#
# Same case and discretization as main.jl. Writes output/convergence.csv with
# columns  CFL, dt, error(2 iterations), error(1 iteration).
#==============================================================================#

using DelimitedFiles

include("scheme.jl")

const OUT_DIR = joinpath(@__DIR__, "output")
mkpath(OUT_DIR)

# Reference: the derived scheme at a dt far below the sweep
CFL_ref = 1/64
println("reference run at CFL $CFL_ref ...")
phi_ref = run_case(CFL_ref, 2; snapshots = 0).phi

CFLs = [1.0, 0.5, 0.25, 0.125, 0.0625]
err2 = Float64[]; err1 = Float64[]
for c in CFLs
    push!(err2, maximum(abs.(run_case(c, 2; snapshots = 0).phi .- phi_ref)))
    push!(err1, maximum(abs.(run_case(c, 1; snapshots = 0).phi .- phi_ref)))
    println("CFL $(c):  dt = $(c*DX)   2 iterations: $(err2[end])   1 iteration: $(err1[end])")
end

# Observed order from a least-squares fit of log(error) against log(dt)
slope(e) = (lx = log.(CFLs .* DX); ly = log.(e); n = length(lx);
            (n*sum(lx .* ly) - sum(lx)*sum(ly)) / (n*sum(lx .^ 2) - sum(lx)^2))
println("fitted order, 2 iterations: ", round(slope(err2), digits=2))
println("fitted order, 1 iteration:  ", round(slope(err1), digits=2))

writedlm(joinpath(OUT_DIR, "convergence.csv"), [CFLs CFLs .* DX err2 err1], ',')
println("output/convergence.csv written (CFL, dt, err 2 iter, err 1 iter)")
