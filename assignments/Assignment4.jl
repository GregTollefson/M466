
using Symbolics
using Printf

# --------------------------------------------------
# Polynomial definitions
# --------------------------------------------------

function p(x)
    prod(x - k for k in 1:10)
end

function q(x, epsilon)
    p(x) + epsilon*x^9
end


# --------------------------------------------------
# Secant method
# --------------------------------------------------

function secant(f, x0, x1; tol=1e-12, maxiter=100)

    for n in 1:maxiter
        f0 = f(x0)
        f1 = f(x1)

        denominator = f1 - f0

        if denominator == 0
            error("Secant method: zero denominator")
        end

        x2 = x1 - f1*(x1-x0)/denominator

        if abs(x2-x1) < tol
            return x2
        end

        x0, x1 = x1, x2
    end

    error("Secant method did not converge")
end


# --------------------------------------------------
# Problem 1(iv): Positive perturbations
# --------------------------------------------------

function problem1_iv()

    alpha = -55.0
    x_previous = 10.0

    results = []

    @printf("%14s %18s %20s\n",
            "epsilon", "beta", "x_beta")

    for k in -24:-13

        epsilon = 2.0^k
        beta = alpha + epsilon

        # Track the root originating at x = 10
        f(x) = q(x, epsilon)

        x_beta = secant(f, x_previous,
                        x_previous - 1e-5)

        @printf("%14.8e %18.12f %20.12f\n",
                epsilon, beta, x_beta)

        push!(results,
              (epsilon=epsilon, beta=beta, x_beta=x_beta))

        x_previous = x_beta
    end

    return results
end


# --------------------------------------------------
# Problem 1(v): Negative perturbations
# --------------------------------------------------

function problem1_v()

    alpha = -55.0
    x_previous = 10.0

    results = []

    @printf("%14s %18s %20s\n",
            "epsilon", "beta", "x_beta")

    for k in -24:-13

        epsilon = -2.0^k
        beta = alpha + epsilon

        # Track the root originating at x = 10
        f(x) = q(x, epsilon)

        x_beta = secant(f, x_previous,
                        x_previous + 1e-5)

        @printf("%14.8e %18.12f %20.12f\n",
                epsilon, beta, x_beta)

        push!(results,
              (epsilon=epsilon, beta=beta, x_beta=x_beta))

        x_previous = x_beta
    end

    return results
end


# --------------------------------------------------
# Problem 2(ii): Third order/Chebyshev Method
# --------------------------------------------------

function problem2_ii()

    setprecision(BigFloat, 16384) do

        # Function and derivatives
        f(x)   = 5*x*exp(-x) - 1
        df(x)  = 5*exp(-x)*(1-x)
        d2f(x) = 5*exp(-x)*(x-2)

        # Chebyshev iteration
        function chebyshev(x)
            fx = f(x)
            dfx = df(x)
            d2fx = d2f(x)

            x - fx/dfx - (fx^2*d2fx)/(2*dfx^3)
        end

        # Initial approximation x_0
        xs = BigFloat[BigFloat(2)]

        # Compute x_1 through x_12
        for n in 1:12
            push!(xs, chebyshev(xs[end]))
        end

        println("\nProblem 2(ii): 3rd order Newtons Method.")
        println("Precision: 16384 bits")
        println("Initial guess: x0 = 2")

        println("\nConvergence estimates:")
        @printf("%4s %16s\n", "n", "p_n")
        println("-------------------------")

        for n in 2:2:10

            d1 = abs(xs[n+1] - xs[n+2])
            d2 = abs(xs[n+2] - xs[n+3])

            if d1 == 0 || d2 == 0
                @printf("%4d %16s\n",
                        n, "precision limit")
                continue
            end

            pn = log(d2)/log(d1)

            @printf("%4d %16.12f\n",
                    n, Float64(pn))
        end

        println("\nApproximate root:")
        @printf("%.30f\n", Float64(xs[end]))

        println("\nFinal residual:")
        @printf("%.6e\n",
                Float64(abs(f(xs[end]))))

        return xs
    end
end


# --------------------------------------------------
# Execute assignments
# --------------------------------------------------

@printf("Running Problem 1(iv)...\n")
results = problem1_iv()

@printf("\nRunning Problem 1(v)...\n")
results_negative = problem1_v()

@printf("\nRunning Problem 2(ii)...\n")
results_2 = problem2_ii()
