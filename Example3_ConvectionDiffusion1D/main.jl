#==============================================================================#
# 1D steady convection-diffusion with Dirichlet BCs.
# a*u_x - nu*u_xx = f,  u(0)=u(1)=0
# CG, linear elements, plain Galerkin (no stabilization)
# Exact solution: u = sin(pi*x)
#==============================================================================#

using Finch

initFinch("convdiff1d");

# results are written to output/, a folder kept in the repo
useLog("convdiff1dlog", dir="output", level=3)

domain(1)
functionSpace(order=1)

mesh(LINEMESH, elsperdim=200, bids=2)

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, 0)
boundary(u, 2, DIRICHLET, 0)

# Physical parameters
coefficient("nu", 0.01)      # diffusivity
coefficient("a", 1.0)        # convection velocity

# Manufactured source: f = a*pi*cos(pi*x) + nu*pi^2*sin(pi*x)
coefficient("f", "1.0*pi*cos(pi*x) + 0.01*pi*pi*sin(pi*x)")

# Weak form: -nu*u_x*v_x - a*u_x*v + f*v = 0
weakForm(u, "-nu*grad(u)*grad(v) - a*grad(u)*v + f*v")

solve(u)

finalizeFinch()

# Error check against exact solution sin(pi*x)
allerr = zeros(size(Finch.finch_state.grid_data.allnodes, 2))

for i = 1:size(Finch.finch_state.grid_data.allnodes, 2)
    local x = Finch.finch_state.grid_data.allnodes[1, i]
    exact = sin(pi*x)
    allerr[i] = abs(u.values[i] - exact)
end

maxerr = maximum(abs, allerr)
println("max error = " * string(maxerr))

using DelimitedFiles

x = Finch.finch_state.grid_data.allnodes[1, :]
uvals = u.values[:]

data = [x uvals]
writedlm("output/solution.csv", data, ',')