using Finch

initFinch("heat2d");

# results are written to output/, a folder kept in the repo
useLog("heat2dlog", dir="output", level=3)

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

boundary(u, 1, DIRICHLET, "0")  # boundary condition for BID 1 is Dirichlet with value 0
boundary(u, 2, DIRICHLET, "0")  # boundary condition for BID 2 is Dirichlet with value 0
boundary(u, 3, DIRICHLET, "0")  # boundary condition for BID 3 is Dirichlet with value 0
boundary(u, 4, DIRICHLET, "1")  # boundary condition for BID 4 is Dirichlet with value 0

coefficient("f", "0")

weakForm(u, "dot(grad(u),grad(v)) - f*v")


solve(u);

finalizeFinch()

using DelimitedFiles

# Get coordinates and solution values
x = Finch.finch_state.grid_data.allnodes[1, :]   # x-coordinates
y = Finch.finch_state.grid_data.allnodes[2, :]   # y-coordinates
uvals = u.values[:]                              # computed solution

# Combine into a 3-column array: x, y, u
data = [x y uvals]

# Write to CSV
writedlm("output/solution.csv", data, ',')






