using DelimitedFiles, Plots, LinearAlgebra

data=readdlm("DE229EC9.dat")
x=data[:,1]
y=data[:,2]
display(scatter(x,y,label="data"))

phi=[x->1,x->x,x->x^2,x->x^3]
V=[phi[j](x[i]) for i=1:length(x), j=1:4]

cfit=V\y
yfit=V*cfit
display(plot!(x,yfit,label="fit"))

ctrue = [2.0, -3.0, 0.5, 2.0]

error = norm(cfit-ctrue)/norm(ctrue)

println("Fit Error = ", error)

ytrue=V*ctrue

display(plot!(x,ytrue,label="true"))
