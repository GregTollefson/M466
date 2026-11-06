using Printf

n = 3

# Map [-1, 1] onto [0, 3]
ell(x) = 1.5 * (x + 1)

@printf("%3s %12s %20s %20s\n",
        "k", "Angle", "Chebyshev node", "Mapped node")

println("-"^62)

for k in 0:n
    theta = (2*k + 1) * pi / (2*(n + 1))

    t = cos(theta)
    x = ell(t)

    @printf("%3d %9d*pi/%-2d %20.6f %20.6f\n",
            k, 2*k + 1, 2*(n + 1), t, x)
end


using Printf
using LinearAlgebra

# Part (i): Affine mapping [-1,1] to [0,3]
ell(x) = 1.5 * (x + 1)

# Function to interpolate
f(x) = sin(pi * x / 6)

# Part (ii): Calculate Chebyshev nodes
n = 3

nodes = [ell(cos((2*k + 1)*pi / (2*(n+1))))
         for k in 0:n]

values = f.(nodes)

println("Chebyshev nodes and function values")
@printf("%3s %14s %14s\n", "k", "x_k", "f(x_k)")
println("-"^35)

for k in 0:n
    @printf("%3d %14.8f %14.8f\n",
            k, nodes[k+1], values[k+1])
end

# Part (iii): Construct cubic interpolating polynomial
# q(x) = c[1] + c[2]*x + c[3]*x^2 + c[4]*x^3

V = [nodes[i]^(j-1) for i in 1:4, j in 1:4]

c = V \ values

println("\nPolynomial coefficients:")
for j in 1:4
    @printf("c[%d] = %.10f\n", j, c[j])
end

# Evaluate q at any point x
q(x) = sum(c[j] * x^(j-1) for j in 1:4)

# Evaluate derivative q'(x)
qprime(x) = sum((j-1)*c[j]*x^(j-2) for j in 2:4)

# Integrate q from 0 to 3
integral_q = sum(c[j]*3^j/j for j in 1:4)

println("\nApproximations:")
@printf("q(0.5)       = %.8f\n", q(0.5))
@printf("q'(1)        = %.8f\n", qprime(1.0))
@printf("Integral 0-3 = %.8f\n", integral_q)

# Exact values for comparison
exact_value = f(0.5)
exact_derivative = (pi/6)*cos(pi/6)
exact_integral = 6/pi

println("\nComparison with exact values:")
@printf("%-14s %14s %14s %14s\n",
        "Quantity", "Approximation", "Exact", "Abs Error")

@printf("%-14s %14.8f %14.8f %14.8f\n",
        "f(0.5)", q(0.5), exact_value,
        abs(q(0.5)-exact_value))

@printf("%-14s %14.8f %14.8f %14.8f\n",
        "f'(1)", qprime(1.0), exact_derivative,
        abs(qprime(1.0)-exact_derivative))

@printf("%-14s %14.8f %14.8f %14.8f\n",
        "Integral", integral_q, exact_integral,
        abs(integral_q-exact_integral))



# ============================================================
# Question 2(v): Plot p(x), q(x), and f(x) on [0,3]
# ============================================================

using Plots

# Equally spaced interpolation nodes
xp = [0.0, 1.0, 2.0, 3.0]
yp = f.(xp)

# Lagrange interpolation for equally spaced nodes
function lagrange(x, nodes, values)
    result = 0.0

    for i in eachindex(nodes)
        L = 1.0

        for j in eachindex(nodes)
            if i != j
                L *= (x - nodes[j]) / (nodes[i] - nodes[j])
            end
        end

        result += values[i] * L
    end

    return result
end

# Polynomial p(x) from Question 1
p(x) = lagrange(x, xp, yp)

# Plot over [0,3]
xplot = range(0, 3, length=500)

plt = plot(xplot, f.(xplot),
    label="f(x) = sin(pi*x/6)",
    linewidth=3,
    color=:black,
    xlabel="x",
    ylabel="y",
    title="Cubic Interpolation Comparison",
    legend=:bottomright)

plot!(plt, xplot, p.(xplot),
    label="p(x): Equally spaced",
    linewidth=2,
    linestyle=:dash,
    color=:blue)

plot!(plt, xplot, q.(xplot),
    label="q(x): Chebyshev",
    linewidth=2,
    linestyle=:dot,
    color=:red)

# Show interpolation nodes
scatter!(plt, xp, yp,
    label="Equally spaced nodes",
    markersize=4,
    color=:blue)

scatter!(plt, nodes, values,
    label="Chebyshev nodes",
    markersize=4,
    color=:red,
    markershape=:diamond)

display(plt)

# Save the plot as a PDF
savefig(plt, "Assignment3_plot.pdf")
