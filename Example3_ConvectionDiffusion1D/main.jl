#==============================================================================#
# 1D steady convection-diffusion with Dirichlet BCs.
# a*u_x - nu*u_xx = f,  u(0)=u(1)=0
# CG, linear elements, plain Galerkin (no stabilization). Finch V2 (branch dev). Run with ./run.sh; the target
# (CPU Julia, or MLIR on the CPU/GPU) is chosen by FINCH_TARGET, see ../common.jl.
# Exact solution: u = sin(pi*x)
#==============================================================================#

include(joinpath(@__DIR__, "..", "common.jl"))   # Finch, FINCH_TARGET switch, helpers

# results are written to output/, a folder kept in the repo; Finch's generated
# project goes to output/finch/
const OUT_DIR  = joinpath(@__DIR__, "output")
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

# A small steady problem: direct solve (LU; the convection term makes it non-symmetric)
ctx = init_finch("convdiff1d"; finch_solver_kwargs(:lu)...,
                 target = finch_target("convdiff1d", PROJ_DIR))

domain(1)
mesh(MESH_LINE; elsperdim = 200)
function_space(order = 1)

u = variable("u")
test_symbol("v")

boundary(u, 1, BC_DIRICHLET, 0.0)
boundary(u, 2, BC_DIRICHLET, 0.0)

# Physical parameters
coefficient("nu", 0.01)      # diffusivity
coefficient("a", 1.0)        # convection velocity

# Manufactured source: f = a*pi*cos(pi*x) + nu*pi^2*sin(pi*x)
coefficient("f", (x, y, z, t) -> 1.0*pi*cos(pi*x) + 0.01*pi*pi*sin(pi*x))

# Weak form: -nu*u_x*v_x - a*u_x*v + f*v = 0
weak_form(u, "-nu*grad(u)*grad(v) - a*grad(u)*v + f*v")

solve(u)
finch_run!(ctx, u)

# Error check against exact solution sin(pi*x)
x = ctx.grid.nodes[1, :]
uvals = u.values[1, :]
maxerr = maximum(abs.(uvals .- sin.(pi .* x)))
println("max error = " * string(maxerr))

data = [x uvals]
writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')
