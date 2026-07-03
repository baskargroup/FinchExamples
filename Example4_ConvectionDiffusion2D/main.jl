#==============================================================================#
# 2D steady convection-diffusion, Dirichlet BCs, plain Galerkin (no SUPG)
# b.grad(u) - nu*div(grad(u)) = f  on [0,1]^2,  u = 0 on boundary
# b = (1,1), nu = 0.1
# Exact solution: u = sin(pi*x)*sin(pi*y)
#==============================================================================#

using Finch

initFinch("convdiff2d");

useLog("convdiff2dlog", level=3)

# Set up the configuration
domain(2)              # dimension
functionSpace(order=2) # basis polynomial order

linAlgOptions(iterative=true)

mesh(QUADMESH, elsperdim=20) # 20x20 elements

addBoundaryID(1, "x==0")
addBoundaryID(2, "x==1")
addBoundaryID(3, "y==0")
addBoundaryID(4, "y==1")

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, "0")
boundary(u, 2, DIRICHLET, "0")
boundary(u, 3, DIRICHLET, "0")
boundary(u, 4, DIRICHLET, "0")

# Physical parameters
coefficient("nu", 0.1)                          # diffusivity
b = coefficient("b", ["1.0", "1.0"], type=VECTOR)  # convection velocity (bx, by)

# Manufactured source for u = sin(pi*x)*sin(pi*y):
# f = bx*pi*cos(pi*x)*sin(pi*y) + by*pi*sin(pi*x)*cos(pi*y) + 2*nu*pi^2*sin(pi*x)*sin(pi*y)
coefficient("f", "1.0*pi*cos(pi*x)*sin(pi*y) + 1.0*pi*sin(pi*x)*cos(pi*y) + 2*0.1*pi*pi*sin(pi*x)*sin(pi*y)")

# Weak form: -nu*grad(u).grad(v) - (b.grad(u))*v + f*v = 0
weakForm(u, "-nu*dot(grad(u),grad(v)) - dot(b,grad(u))*v + f*v")

solve(u);

finalizeFinch()

# Error check against exact solution
maxerr = 0.0
for i = 1:size(Finch.finch_state.grid_data.allnodes, 2)
    local xn = Finch.finch_state.grid_data.allnodes[1, i]
    local yn = Finch.finch_state.grid_data.allnodes[2, i]
    exact = sin(pi*xn)*sin(pi*yn)
    global maxerr = max(maxerr, abs(u.values[i] - exact))
end
println("max error = " * string(maxerr))

using DelimitedFiles

# Get coordinates and solution values
x = Finch.finch_state.grid_data.allnodes[1, :]
y = Finch.finch_state.grid_data.allnodes[2, :]
uvals = u.values[:]

data = [x y uvals]
writedlm("solution.csv", data, ',')