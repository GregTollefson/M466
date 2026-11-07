
# residual.jl
# Compute residual for the discrete Poisson equation

using LinearAlgebra

function compute_residual(u, f, h)
    Au = apply_operator(u, h)
    return f - Au
end

function relative_residual(u, f, h)
    r = compute_residual(u, f, h)
    return norm(r) / norm(f)
end
