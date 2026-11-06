
using Symbolics

@variables x

# Rodrigues formula
function laguerre(n)
    expr = x^n * exp(-x)

    for k in 1:n
        expr = expand_derivatives(Differential(x)(expr))
    end

    return simplify(exp(x) * expr)
end

# Calculate p3, p4, p5
for n in 3:5
    println("p_$n(x) = ", laguerre(n))
end
