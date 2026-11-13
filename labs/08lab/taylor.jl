# taylor.jl -- Find the Taylor polynomial coefficients of degree n

using Symbolics

function taylor(p, n)
    @variables t

    D(g) = expand_derivatives(Differential(t)(g))

    dp = p(t)
    a = zeros(n + 1)

    for k = 1:n+1
        a[k] = Symbolics.value(substitute(dp, t => 0))
        dp = D(dp / k)
    end

    return a
end