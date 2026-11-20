
using LinearAlgebra

# Construct the n x n Jacobi matrix
function makeJn(n)
    Jn = zeros(n, n)

    for j = 1:n-1
        z = j / sqrt(4*j^2 - 1)

        Jn[j+1, j] = z
        Jn[j, j+1] = z
    end

    return Jn
end

# Test with a 4 x 4 matrix
J4 = makeJn(4)
display(J4)


# Compute eigenvalues and eigenvectors
xn, U = eigen(J4)

# Eigenvalues are the quadrature nodes
println("Nodes:")
display(xn)

# Calculate weights from the first row of U
wn = 2 * U[1, :].^2

println("Weights:")
display(wn)

