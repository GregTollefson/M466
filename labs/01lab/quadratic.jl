function xplus(a,b,c)
    return (-b + sqrt(b^2 - 4*a*c)) / (2*a)
end
function xminus(a,b,c)
    return (-b - sqrt(b^2 - 4*a*c)) / (2*a)
end
function xplusnew(a,b,c)
    return 2c/(-b - sqrt(b^2 - 4*a*c))
end
function xminusnew(a,b,c)
    return 2c/(-b + sqrt(b^2 - 4*a*c))
end
println("xplus = ", xplus(1,1,-5))
println("xminus = ", xminus(1,1,-5))
println("xplusnew = ", xplusnew(1,1,-5))
println("xminusnew = ", xminusnew(1,1,-5))