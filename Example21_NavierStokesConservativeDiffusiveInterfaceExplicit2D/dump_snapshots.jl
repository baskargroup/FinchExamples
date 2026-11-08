# Write the snapshots of a run (also one that blew up before main.jl wrote its
# csv files) as the same alpha/ux/uy/p_timeseries.csv and snapshot_times.csv.
# Usage: CDI_NEL=128 CDI_LY=0.75 julia +1.12 --project=$FINCH_V2_DIR dump_snapshots.jl output_caseD_halfdt
using Finch
using DelimitedFiles
out = joinpath(@__DIR__, ARGS[1])
NEL = parse(Int, get(ENV, "CDI_NEL", "128")); LY = parse(Float64, get(ENV, "CDI_LY", "0.75"))
ctx = init_finch("dump"; verbose = false, target = JuliaCPUTarget(out_dir = joinpath(out, "finch"), func_name = "dump"))
domain(2)
m = mesh(MESH_QUAD; elsperdim = (NEL, round(Int, LY * NEL)))
m.nodes[2, :] .*= LY
function_space(order = 1)
g = ctx.grid === nothing ? m : ctx.grid      # order 1: the FE nodes are the mesh nodes
x = g.nodes[1, :]; y = g.nodes[2, :]; N = length(x)
results = joinpath(out, "finch", "results")
_, snaps = Finch.read_results_manifest(results)
vals  = [Finch.load_snapshot_values(results, s) for s in snaps]
times = [Float64(s["t"]) for s in snaps]
A  = reduce(hcat, [v["alpha"] for v in vals])
UX = reduce(hcat, [v["u"][1:N] for v in vals])
UY = reduce(hcat, [v["u"][N+1:2N] for v in vals])
P  = reduce(hcat, [v["p"] for v in vals])
writedlm(joinpath(out, "alpha_timeseries.csv"), [x y A], ',')
writedlm(joinpath(out, "ux_timeseries.csv"), [x y UX], ',')
writedlm(joinpath(out, "uy_timeseries.csv"), [x y UY], ',')
writedlm(joinpath(out, "p_timeseries.csv"), [x y P], ',')
writedlm(joinpath(out, "snapshot_times.csv"), times, ',')
println(length(times), " snapshots written, N = ", N, ", last t = ", times[end])
