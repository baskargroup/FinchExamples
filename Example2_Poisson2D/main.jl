#==============================================================================#
# 2D Laplace problem on the unit square: u = 1 on the top side, 0 on the others.
# CG, quadratic elements. Finch V2, Julia target.
#==============================================================================#

using Finch
using DelimitedFiles

# results are written to output/, a folder kept in the repo; Finch's generated
# project goes to output/finch/
const OUT_DIR  = joinpath(@__DIR__, "output")
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

ctx = init_finch("heat2d"; solver = :lu,
                 target = JuliaCPUTarget(out_dir = PROJ_DIR, func_name = "heat2d"))

# Set up the configuration
domain(2)                          # dimension
mesh(MESH_QUAD; elsperdim = 20)    # 20x20 elements
function_space(order = 2)          # basis polynomial order

# Boundary ids: 1 = x==0, 2 = x==1, 3 = y==0, 4 = y==1

u = variable("u")
test_symbol("v")

boundary(u, 1, BC_DIRICHLET, 0.0)  # boundary condition for BID 1 is Dirichlet with value 0
boundary(u, 2, BC_DIRICHLET, 0.0)  # boundary condition for BID 2 is Dirichlet with value 0
boundary(u, 3, BC_DIRICHLET, 0.0)  # boundary condition for BID 3 is Dirichlet with value 0
boundary(u, 4, BC_DIRICHLET, 1.0)  # boundary condition for BID 4 is Dirichlet with value 1

coefficient("f", 0.0)

weak_form(u, "dot(grad(u),grad(v)) - f*v")

solve(u)
Finch.run_emitted_julia!(ctx.target, u)

# Get coordinates and solution values
x = ctx.grid.nodes[1, :]   # x-coordinates
y = ctx.grid.nodes[2, :]   # y-coordinates
uvals = u.values[1, :]     # computed solution

# Combine into a 3-column array: x, y, u
data = [x y uvals]

# Write to CSV
writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')
