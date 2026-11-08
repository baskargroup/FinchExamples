#==============================================================================#
# 1D Poisson problem with Dirichlet boundary conditions.
# CG, Linear element. Finch V2 (branch dev). Run with ./run.sh; the target
# (CPU Julia, or MLIR on the CPU/GPU) is chosen by FINCH_TARGET, see ../common.jl.
#==============================================================================#

include(joinpath(@__DIR__, "..", "common.jl"))   # Finch, FINCH_TARGET switch, helpers

# results are written to output/, a folder kept in the repo; Finch's generated
# project goes to output/finch/
const OUT_DIR  = joinpath(@__DIR__, "output")
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

# Initialize 1D Poisson problem setup
ctx = init_finch("poisson1d"; finch_solver_kwargs(:lu)...,
                 target = finch_target("poisson1d", PROJ_DIR))

domain(1)                       # Define problem as 1-dimensional
mesh(MESH_LINE; elsperdim = 200)   # Create uniform 1D mesh with 200 elements
function_space(order = 1)       # polynomial order

u = variable("u")               # Define scalar unknown field u
test_symbol("v")                # Set test function symbol to v

# Apply Dirichlet boundary condition
boundary(u, 1, BC_DIRICHLET, 0.0)
boundary(u, 2, BC_DIRICHLET, 0.0)

# Define source term f(x)
coefficient("f", (x, y, z, t) -> -100*pi*pi*sin(10*pi*x)*sin(pi*x) - pi*pi*sin(10*pi*x)*sin(pi*x) +
                                 20*pi*pi*cos(10*pi*x)*cos(pi*x))

# Define weak form of PDE: ∫ (∇u·∇v + f*v) dx = 0
weak_form(u, "-grad(u)*grad(v) - f*v")

# Solve the PDE
solve(u)
finch_run!(ctx, u)

# Exact solution: sin(10*pi*x)*sin(pi*x)
x = ctx.grid.nodes[1, :]                     # x-coordinates
_, fields = read_snapshots(PROJ_DIR)        # the solution as both targets write it
uvals = final_values(fields, "u")                       # computed solution
allerr = abs.(uvals .- sin.(10*pi .* x) .* sin.(pi .* x))   # absolute error at every node

maxerr = maximum(allerr)                     # Compute maximum absolute error
println("max error = " * string(maxerr))     # Print it

data = [x uvals]                             # two-column array
writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')   # write to CSV
