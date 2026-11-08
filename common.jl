#==============================================================================#
# Shared pieces for every example: the FINCH_TARGET switch, running the emitted
# project, and the post-processing blocks the CDI examples share.
#
# Include from an example's main.jl:   include(joinpath(@__DIR__, "..", "common.jl"))
#
# FINCH_TARGET (environment variable):
#   julia      (default) Finch's Julia target: the solve runs in this process (CPU, the reference).
#   mlir_cpu   Finch's MLIR/IREE target on the CPU: the script emits a project, then runs its
#              build.sh and run.sh (needs IREE: `setup_iree()` once in the Finch checkout, Linux).
#   mlir_cuda  the same on an NVIDIA GPU (CUDA toolkit and nvidia-smi on the machine, or FINCH_CUDA_ARCH).
# The physics (mesh, variables, weak forms, schedule) is written once; only the target, the
# solver settings the MLIR target requires, and the way the run is launched differ.
#==============================================================================#

using Finch
using DelimitedFiles

const FINCH_TARGETS = ("julia", "mlir_cpu", "mlir_cuda")
const FINCH_TARGET  = get(ENV, "FINCH_TARGET", "julia")
FINCH_TARGET in FINCH_TARGETS ||
    error("FINCH_TARGET must be julia, mlir_cpu or mlir_cuda, got \"$FINCH_TARGET\"")

"""
    finch_target(func_name, proj_dir)

The Finch code target for `FINCH_TARGET`, emitting into `proj_dir`.
"""
function finch_target(func_name::String, proj_dir::String)
    FINCH_TARGET == "julia" && return JuliaCPUTarget(out_dir = proj_dir, func_name = func_name)
    device = FINCH_TARGET == "mlir_cuda" ? :cuda : :cpu
    return MLIRTarget(out_dir = proj_dir, func_name = func_name, device = device)
end

"""
    finch_solver_kwargs(julia_solver; tol = 1e-9, max_iter = 500)

Keywords for `init_finch`: the example's own solver on the Julia target; BiCGSTAB with Jacobi
on the MLIR target (its matrix-free solvers are CG and BiCGSTAB only).
"""
function finch_solver_kwargs(julia_solver::Symbol; tol::Float64 = 1e-9, max_iter::Int = 500)
    FINCH_TARGET == "julia" &&
        return (solver = julia_solver, solver_tol = tol, solver_max_iter = max_iter)
    return (solver = :bicgstab, preconditioner = :jacobi, solver_tol = tol, solver_max_iter = max_iter)
end

"""
    finch_run!(ctx, vars...)

Run the project `solve()` emitted: in this process on the Julia target (fills `var.values`);
through the project's build.sh and run.sh on the MLIR target (results land in `results/`).
"""
function finch_run!(ctx, vars...)
    if FINCH_TARGET == "julia"
        Finch.run_emitted_julia!(ctx.target, vars...)
    else
        proj_dir = ctx.target.out_dir
        run(Cmd(`bash build.sh`; dir = proj_dir))
        run(Cmd(`bash run.sh`; dir = proj_dir))
    end
    return nothing
end

# ---- results written by both targets ------------------------------------------

"""
    read_snapshots(proj_dir) -> (times, fields)

All snapshots of the run: `times[k]` and `fields[name][:, k]` (nodes x snapshots; a vector field
is component-blocked, see `final_vector`).
"""
function read_snapshots(proj_dir::String)
    results = joinpath(proj_dir, "results")
    _, snaps = Finch.read_results_manifest(results)
    vals  = [Finch.load_snapshot_values(results, s) for s in snaps]
    times = [Float64(s["t"]) for s in snaps]
    names = isempty(vals) ? String[] : collect(keys(vals[1]))
    fields = Dict{String, Matrix{Float64}}(n => reduce(hcat, [v[n] for v in vals]) for n in names)
    return times, fields
end

final_values(fields, name::String) = fields[name][:, end]

"""
    check_final_time(times, T)

Error unless the last snapshot is the final time `T` (both targets write the final step as a
snapshot; the post-processing of every example reads the final state from it).
"""
function check_final_time(times, T)
    abs(times[end] - T) <= 1e-9 * T ||
        error("the last snapshot is at t = $(times[end]), not the final time T = $T")
    return nothing
end
final_vector(fields, name::String, N::Int) = (fields[name][1:N, end], fields[name][N+1:2N, end])

write_timeseries(out_dir::String, name::String, coords::AbstractMatrix, A::AbstractMatrix) =
    writedlm(joinpath(out_dir, name * "_timeseries.csv"), [coords A], ',')
write_snapshot_times(out_dir::String, times) =
    writedlm(joinpath(out_dir, "snapshot_times.csv"), times, ',')

"""
    lumped_weights(x, y, dx)

Nodal quadrature weights of a uniform quad mesh on [0,1]^2 with spacing `dx`: dx^2 inside,
halved on an edge, quartered at a corner (the L1 norm of the explicit and semi-implicit examples).
"""
function lumped_weights(x, y, dx)
    w = fill(dx^2, length(x))
    for i in eachindex(x)
        onx = (x[i] < 1e-12 || x[i] > 1 - 1e-12)
        ony = (y[i] < 1e-12 || y[i] > 1 - 1e-12)
        onx && (w[i] /= 2)
        ony && (w[i] /= 2)
    end
    return w
end

"""
    mass_change(results_dir, col = 3) -> (m0, m1, relative)

First and last value of a monitor column of `monitors.csv` (column 3 is the first monitor).
"""
function mass_change(results_dir::String, col::Int = 3)
    mon = readdlm(joinpath(results_dir, "monitors.csv"), ',', skipstart = 1)
    m0, m1 = mon[1, col], mon[end, col]
    return m0, m1, (m1 - m0) / m0
end

# ---- the CDI normal n = grad(phi)/|grad(phi)| contracted with grad(v) -----------
# nodal (recovered) normal through the projected gradient g (or gx in 1D), and the
# element-wise normal straight from grad(phi); delta regularises |grad phi| = 0.
const NDOTGV_NODAL_1D   = "gx*pow(gx*gx + delta, -0.5)"
const NDOTGV_NODAL_2D   = "(g[1]*grad(v)[1] + g[2]*grad(v)[2])*pow(g[1]*g[1] + g[2]*g[2] + delta, -0.5)"
const NDOTGV_ELEMENT_2D = "dot(grad(phi), grad(v))*pow(dot(grad(phi), grad(phi)) + delta, -0.5)"
