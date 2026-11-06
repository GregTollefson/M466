
setprecision(4096)

#f(x) = 5*x*exp(-x) - 1
#df(x) = 5*(1-x)*exp(-x)
#x0 = big"2.0"

f(x) = cos(x) - 2*x
df(x) = -sin(x) - 2

x0 = big"0.5" 

g(x) = x - f(x)/df(x)

xn = x0
xn1 = g(xn)
en = abs(xn - xn1)

for n = 1:9
    global xn, xn1, en

    xn = xn1
    xn1 = g(xn)

    en0 = en
    en = abs(xn - xn1)

    println("n=$(n-1) log(e(n+1))/log(e(n))=",
            Float64(log(en)/log(en0)))
end
