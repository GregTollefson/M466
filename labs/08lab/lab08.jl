using SpecialMatrices, Polynomials, LinearAlgebra

include("taylor.jl")

# Define the polynomial
p(x) = x^5 - 9*x^4 - 4*x^3 - 6*x^2 + 4

# Compute polynomial coefficients
a = taylor(p, 5)

# Construct the companion matrix
Cp = Companion(Polynomial(a))

# Compute the eigenvalues (roots)
lambda = eigvals(Cp)

# Verify the roots
for x in lambda
    println("|p($x)| = ", abs(p(x)))
end


# Construct the taylor polynomial for differentiation
poly = Polynomial(a)

# Compute derivative polynomial
dp = derivative(poly)

# Calculate residual and estimated root error
for x in lambda
    residual = abs(p(x))
    error_est = abs(p(x) / dp(x))

    println("Root: ", x)
    println("  Residual:        ", residual)
    println("  Estimated error: ", error_est)
    println()
end