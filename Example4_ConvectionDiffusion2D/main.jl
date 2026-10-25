#==============================================================================#
# 2D steady convection-diffusion, Dirichlet BCs, plain Galerkin (no SUPG)
# b.grad(u) - nu*div(grad(u)) = f  on [0,1]^2,  u = 0 on boundary
# b = (1,1), nu = 0.1
# Exact solution: u = sin(pi*x)*sin(pi*y)
# Finch V2, Julia target.
#==============================================================================#

using Finch
using DelimitedFiles

# results are written to output/, a folder kept in the repo; Finch's generated
# project goes to output/finch/
const OUT_DIR  = joinpath(@__DIR__, "output")
const PROJ_DIR = joinpath(OUT_DIR, "finch")
mkpath(PROJ_DIR)

# A small steady problem: direct solve (LU; the convection term makes it non-symmetric)
ctx = init_finch("convdiff2d"; solver = :lu,
                 target = JuliaCPUTarget(out_dir = PROJ_DIR, func_name = "convdiff2d"))

# Set up the configuration
domain(2)                          # dimension
mesh(MESH_QUAD; elsperdim = 20)    # 20x20 elements
function_space(order = 2)          # basis polynomial order

# Boundary ids: 1 = x==0, 2 = x==1, 3 = y==0, 4 = y==1

u = variable("u")
test_symbol("v")

for bid in 1:4
    boundary(u, bid, BC_DIRICHLET, 0.0)
end

# Physical parameters
coefficient("nu", 0.1)                                          # diffusivity
coefficient("b", [(x, y, z, t) -> 1.0, (x, y, z, t) -> 1.0])    # convection velocity (bx, by)

# Manufactured source for u = sin(pi*x)*sin(pi*y):
# f = bx*pi*cos(pi*x)*sin(pi*y) + by*pi*sin(pi*x)*cos(pi*y) + 2*nu*pi^2*sin(pi*x)*sin(pi*y)
coefficient("f", (x, y, z, t) -> 1.0*pi*cos(pi*x)*sin(pi*y) + 1.0*pi*sin(pi*x)*cos(pi*y) +
                                 2*0.1*pi*pi*sin(pi*x)*sin(pi*y))

# Weak form: -nu*grad(u).grad(v) - (b.grad(u))*v + f*v = 0
weak_form(u, "-nu*dot(grad(u),grad(v)) - dot(b,grad(u))*v + f*v")

solve(u)
Finch.run_emitted_julia!(ctx.target, u)

# Error check against exact solution
x = ctx.grid.nodes[1, :]
y = ctx.grid.nodes[2, :]
uvals = u.values[1, :]
maxerr = maximum(abs.(uvals .- sin.(pi .* x) .* sin.(pi .* y)))
println("max error = " * string(maxerr))

data = [x y uvals]
writedlm(joinpath(OUT_DIR, "solution.csv"), data, ',')
