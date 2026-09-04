using Plots

# Test Function
#f(x) =5*x*exp(-x)-1
#x0 = 2

# Assignment Function
f(x) = cos(x) - x/2
x0 = 1

xs=0:0.1:5
ys=f.(xs)

p = plot(xs,ys)
scatter!([x0],[f(x0)])

savefig(p, "graph02.pdf")

gui(p)
