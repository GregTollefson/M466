# Phase 1: CPU Reference Implementation
# grid.jl - Initialize the computational grid

N = 3                  # Number of interior points per dimension
h = 1.0 / (N + 1)      # Grid spacing

# Solution array, including boundary points
u = zeros(Float64, N + 2, N + 2)

# Right-hand side: f(x,y) = 1
f = ones(Float64, N, N)

println("Interior grid size: ", N, " x ", N)
println("Grid spacing h = ", h)
println("Solution array u:")
display(u)

println("Source array f:")
display(f)