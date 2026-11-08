
using LinearAlgebra
using SparseArrays
using Printf

include("operator.jl")
include("redblack.jl")
include("residual.jl")
include("multigrid.jl")

# Grid setup
N = 15
h = 1.0 / (N + 1)

u = zeros(Float64, N+2, N+2)
f = ones(Float64, N, N)

# Step 1: Solve using multigrid
cycles = 0
res = relative_residual(u, f, h)


while res > 1e-8 && cycles < 50
    vcycle!(u, f, h)
    global cycles += 1
    global res = relative_residual(u, f, h)
end


# Extract the interior solution as a vector
u_mg = vec(u[2:N+1, 2:N+1])

# Step 2: Construct the sparse finite-difference matrix
# Julia column-major ordering: k = i + (j-1)*N

e = ones(N)
T = spdiagm(-1 => -e[1:end-1],
             0 =>  2*e,
             1 => -e[1:end-1])

I_N = sparse(1.0I, N, N)

A = (kron(I_N, T) + kron(T, I_N)) / h^2

# Step 3: Direct solution
b = vec(f)
u_direct = A \ b

# Step 4: Compare
relative_error = norm(u_mg - u_direct) / norm(u_direct)
direct_residual = norm(A*u_direct - b) / norm(b)

@printf("Grid size N             = %d\n", N)
@printf("Number of unknowns      = %d\n", N^2)
@printf("Multigrid V-cycles      = %d\n", cycles)
@printf("Multigrid residual      = %.6e\n", res)
@printf("Direct solver residual  = %.6e\n", direct_residual)
@printf("Relative solution error = %.6e\n", relative_error)
