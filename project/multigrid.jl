
# multigrid.jl
# Multigrid components for the Poisson solver

function restrict_full_weighting(r)
    N = size(r, 1)
    @assert N >= 3 && isodd(N)

    Nc = (N - 1) ÷ 2
    rc = zeros(Float64, Nc, Nc)

    # Pad fine residual with zero boundary values
    rp = zeros(Float64, N + 2, N + 2)
    rp[2:N+1, 2:N+1] .= r

    for J in 1:Nc
        for I in 1:Nc
            i = 2*I + 1
            j = 2*J + 1

            rc[I,J] = (
                4*rp[i,j] +
                2*(rp[i-1,j] + rp[i+1,j] +
                   rp[i,j-1] + rp[i,j+1]) +
                rp[i-1,j-1] + rp[i-1,j+1] +
                rp[i+1,j-1] + rp[i+1,j+1]
            ) / 16
        end
    end

    return rc
end


# Bilinear interpolation of coarse-grid correction
# ec contains only interior coarse-grid values.
# Return only interior fine-grid values.




function prolong_bilinear(ec)
    Nc = size(ec, 1)
    Nf = 2*Nc + 1

    # Coarse grid including zero boundaries
    c = zeros(Float64, Nc + 2, Nc + 2)
    c[2:Nc+1, 2:Nc+1] .= ec

    # Fine grid including zero boundaries
    v = zeros(Float64, Nf + 2, Nf + 2)

    # Copy coarse points to fine grid
    for j in 1:Nc+2
        for i in 1:Nc+2
            v[2*i-1, 2*j-1] = c[i,j]
        end
    end

    # Interpolate halfway between coarse points
    # in the first coordinate
    for j in 1:2:Nf+2
        for i in 2:2:Nf+1
            v[i,j] = (v[i-1,j] + v[i+1,j]) / 2
        end
    end

    # Interpolate halfway between coarse points
    # in the second coordinate
    for j in 2:2:Nf+1
        for i in 1:Nf+2
            v[i,j] = (v[i,j-1] + v[i,j+1]) / 2
        end
    end

    return v[2:Nf+1, 2:Nf+1]
end


# One recursive multigrid V-cycle
# u includes boundary points
# f contains interior source values

function vcycle!(u, f, h; nu1=2, nu2=2)

    N = size(f, 1)

    # Coarsest grid: one interior unknown
    if N == 1
        u[2,2] = (h^2 / 4) * f[1,1]
        return u
    end

    @assert N >= 3 && isodd(N)

    # 1. Pre-smoothing
    for k in 1:nu1
        redblack_step!(u, f, h)
    end

    # 2. Compute fine-grid residual
    r = compute_residual(u, f, h)

    # 3. Restrict residual to coarse grid
    rc = restrict_full_weighting(r)

    # 4. Solve coarse-grid error equation
    Nc = size(rc, 1)
    ec = zeros(Float64, Nc+2, Nc+2)

    vcycle!(ec, rc, 2*h; nu1=nu1, nu2=nu2)

    # 5. Prolongate correction
    ef = prolong_bilinear(ec[2:Nc+1, 2:Nc+1])

    # 6. Correct fine-grid approximation
    u[2:N+1, 2:N+1] .+= ef

    # 7. Post-smoothing
    for k in 1:nu2
        redblack_step!(u, f, h)
    end

    return u
end
