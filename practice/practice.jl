println("this is a test")
n = 20
x = rand(1:10,n)
x_mean = sum(x)/n
s_x = sqrt(sum((x .- x_mean).^2)/(n-1)) # sample standard deviation
println("Sample mean: $x_mean")
println("Sample standard deviation: $s_x")