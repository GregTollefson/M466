
# soln04.jl -- evaluate and compute the error in p(2.25)
using Plots
include("lagrange.jl")
x=1.5:0.5:3.5
y=log.(x)
println("p(2.25)=",pp(2.25))
println("E=",abs(pp(2.25)-log(2.25)))
plot([log,pf],1:0.01:4,label=["log" "p"])
scatter!(x,y,label="")
savefig("poly.pdf")
