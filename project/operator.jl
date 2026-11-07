
# Phase 1: CPU Reference Implementation
# operator.jl - Five-point finite-difference operator

function apply_operator(u, h)
    N = size(u, 1) - 2

    Au = zeros(Float64, N, N)

    for j in 2:N+1
        for i in 2:N+1

            Au[i-1, j-1] = (
                4*u[i,j]
                - u[i-1,j]
                - u[i+1,j]
                - u[i,j-1]
                - u[i,j+1]
            ) / h^2

        end
    end

    return Au
end
