
using LinearAlgebra
using Printf

# Construct the Jacobi matrix
function makeJn(n)
    Jn = zeros(n, n)

    for j = 1:n-1
        z = j / sqrt(4*j^2 - 1)
        Jn[j+1, j] = z
        Jn[j, j+1] = z
    end

    return Jn
end

# Function to integrate
f(x) = 1 / (1 + (x - 1)^2)

# Exact integral
Iexact = atan(2)

println("Gaussian Quadrature Results")
@printf("%3s  %20s  %15s\n", "n", "Qn", "En")
println("-"^44)

for n = 1:15
    Jn = makeJn(n)

    # Nodes and weights
    xn, U = eigen(Jn)
    wn = 2 * U[1, :].^2

    # Quadrature approximation
    Qn = f.(xn)' * wn

    # Absolute error
    En = abs(Iexact - Qn)

    @printf("%3d  %20.15f  %15.8e\n", n, Qn, En)
end
