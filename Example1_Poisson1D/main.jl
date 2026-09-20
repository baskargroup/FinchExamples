#==============================================================================#
# 1D Poisson problem with Dirichlet boundary conditions.
# CG, Linear element
#==============================================================================#

using Finch                     # Load Finch package

initFinch("poisson1d");         # Initialize  1D Poisson problem setup


# results are written to output/, a folder kept in the repo
useLog("poisson1dlog", dir="output", level=3) # Enable logging to file

# Set up the configuration

domain(1)                       # Define problem as 1-dimensional
functionSpace(order=1) 		    # polynomial order

mesh(LINEMESH, elsperdim=200, bids=2)   # Create uniform 1D mesh with 200 elements

u = variable("u")               # Define scalar unknown field u

testSymbol("v")                 # Set test function symbol to v


# Apply Dirichlet boundary condition

boundary(u, 1, DIRICHLET, 0)
boundary(u, 2, DIRICHLET, 0)

# Write the weak form 
# Define source term f(x)
coefficient("f", "-100*pi*pi*sin(10*pi*x)*sin(pi*x) - pi*pi*sin(10*pi*x)*sin(pi*x) + 20*pi*pi*cos(10*pi*x)*cos(pi*x)")

# Define weak form of PDE: ∫ (∇u·∇v + f*v) dx = 0
weakForm(u, "-grad(u)*grad(v) - f*v")

# Solve the PDE
solve(u)                   

finalizeFinch()   # Clean up Finch after simulation 


# Exact solution: sin(10*pi*x)*sin(pi*x)
# Prepare array to store error at all nodes
allerr = zeros(size(Finch.finch_state.grid_data.allnodes, 2))

for i = 1:size(Finch.finch_state.grid_data.allnodes, 2)
   local x = Finch.finch_state.grid_data.allnodes[1, i]  # get x-coordinate of node i
    exact = sin(10*pi*x) * sin(pi*x)                # exact solution at x
    allerr[i] = abs(u.values[i] - exact)            # absolute error at node i
end


maxerr = maximum(abs, allerr)               # Compute maximum absolute error
println("max error = " * string(maxerr))    # Print it


using DelimitedFiles

x = Finch.finch_state.grid_data.allnodes[1, :]  # x-coordinates
uvals = u.values[:]                              # computed solution

data = [x uvals]                                # two-column array
writedlm("output/solution.csv", data, ',')             # write to CSV





