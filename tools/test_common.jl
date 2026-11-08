# Tests for common.jl. Run: julia +1.12 --project=$HOME/Packages/Finch-V2 tools/test_common.jl
using Test
const ROOT = normpath(joinpath(@__DIR__, ".."))
delete!(ENV, "FINCH_TARGET")
include(joinpath(ROOT, "common.jl"))

function in_fresh_process(target, code)
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -e $code`
    env = copy(ENV); env["FINCH_TARGET"] = target
    out = IOBuffer(); err = IOBuffer()
    ok = success(pipeline(setenv(cmd, env); stdout = out, stderr = err))
    return ok, String(take!(out)) * String(take!(err))
end
const LOAD = "include(\"$(joinpath(ROOT, "common.jl"))\"); "

@testset "common.jl" begin
    @test FINCH_TARGET == "julia"
    # 3x3 nodes on the unit square, dx = 0.5: corners dx^2/4, edges dx^2/2, centre dx^2 -> total 1 = area
    x = repeat(0.0:0.5:1.0, 3); y = repeat(0.0:0.5:1.0, inner = 3)
    w = lumped_weights(x, y, 0.5)
    @test sum(w) ≈ 1.0
    @test w[5] ≈ 0.25 && w[1] ≈ 0.0625 && w[2] ≈ 0.125
    @test finch_solver_kwargs(:lu).solver == :lu
    @test finch_solver_kwargs(:cg; tol = 1e-12, max_iter = 5000) == (solver = :cg, solver_tol = 1e-12, solver_max_iter = 5000)
    ok, out = in_fresh_process("mlir_cuda", LOAD * "print(finch_solver_kwargs(:lu).solver, ' ', finch_solver_kwargs(:lu).preconditioner)")
    @test ok && occursin("bicgstab jacobi", out)
    ok, out = in_fresh_process("foo", LOAD * "println(\"reached\")")
    @test !ok && occursin("FINCH_TARGET must be julia, mlir_cpu or mlir_cuda", out) && !occursin("reached", out)
    @test check_final_time([0.0, 0.5, 1.0], 1.0) === nothing
    @test_throws ErrorException check_final_time([0.0, 0.5, 0.9], 1.0)
    @test NDOTGV_NODAL_2D == "(g[1]*grad(v)[1] + g[2]*grad(v)[2])*pow(g[1]*g[1] + g[2]*g[2] + delta, -0.5)"
end
