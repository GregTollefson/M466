
# redblack.jl
# One red-black Gauss-Seidel iteration

function redblack_step!(u, f, h)
    N = size(f, 1)

    # Red pass: even parity
    for j in 2:N+1
        for i in 2:N+1
            if iseven(i + j)
                u[i,j] = (
                    u[i-1,j] + u[i+1,j] +
                    u[i,j-1] + u[i,j+1] +
                    h^2 * f[i-1,j-1]
                ) / 4
            end
        end
    end

    # Black pass: odd parity
    for j in 2:N+1
        for i in 2:N+1
            if isodd(i + j)
                u[i,j] = (
                    u[i-1,j] + u[i+1,j] +
                    u[i,j-1] + u[i,j+1] +
                    h^2 * f[i-1,j-1]
                ) / 4
            end
        end
    end

    return u
end


function solve_redblack!(u, f, h;
                         tol=1e-8,
                         maxiter=10000)

    for k in 1:maxiter
        redblack_step!(u, f, h)

        res = relative_residual(u, f, h)

        if res < tol
            println("Converged after ", k, " iterations")
            println("Relative residual = ", res)
            return u
        end
    end

    println("Maximum iterations reached")
    println("Relative residual = ",
            relative_residual(u, f, h))

    return u
end
