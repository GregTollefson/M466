using Random, DelimitedFiles, Plots

Random.seed!(123)

N = 100

# Generate x values
x = collect(range(-2, 2, length=N))

# Construct the 100 × 4 Vandermonde matrix
A = hcat(ones(N), x, x.^2, x.^3)

# True polynomial coefficients
β_true = [2.0, -3.0, 0.5, 2.0]

# Generate normally distributed noise
σ = 1.0
ϵ = σ .* randn(N)

# Generate response vector
y = A * β_true + ϵ

# Combine x and y into a 100 × 2 matrix
data = hcat(x, y)

# Save as a whitespace-delimited text file
writedlm("DE229EC9.dat", data)