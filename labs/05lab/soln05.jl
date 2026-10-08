# soln04.jl -- evaluate and compute the error in p(2.25)

include("lagrange.jl")
x=1.5:0.5:3.5
y=log.(x)
println("p(2.25)=",pp(2.25))
println("E=",abs(pp(2.25)-log(2.25)))
