# Write the snapshots of a run (also one that blew up before main.jl wrote its
# csv files) as the same alpha/ux/uy/p_timeseries.csv and snapshot_times.csv.
# Usage: CDI_NEL=128 CDI_LY=0.75 julia +1.12 --project=$FINCH_V2_DIR dump_snapshots.jl output_caseD_halfdt
include(joinpath(@__DIR__, "..", "common.jl"))
out = joinpath(@__DIR__, ARGS[1])
NEL = parse(Int, get(ENV, "CDI_NEL", "128")); LY = parse(Float64, get(ENV, "CDI_LY", "0.75"))
ctx = init_finch("dump"; verbose = false, target = JuliaCPUTarget(out_dir = joinpath(out, "finch"), func_name = "dump"))
domain(2)
m = mesh(MESH_QUAD; elsperdim = (NEL, round(Int, LY * NEL)))
m.nodes[2, :] .*= LY
function_space(order = 1)
g = ctx.grid === nothing ? m : ctx.grid      # order 1: the FE nodes are the mesh nodes
x = g.nodes[1, :]; y = g.nodes[2, :]; N = length(x)
times, fields = read_snapshots(joinpath(out, "finch"))
A  = fields["alpha"]
UX = fields["u"][1:N, :]
UY = fields["u"][N+1:2N, :]
P  = fields["p"]
write_timeseries(out, "alpha", [x y], A)
write_timeseries(out, "ux", [x y], UX)
write_timeseries(out, "uy", [x y], UY)
write_timeseries(out, "p", [x y], P)
write_snapshot_times(out, times)
println(length(times), " snapshots written, N = ", N, ", last t = ", times[end])
