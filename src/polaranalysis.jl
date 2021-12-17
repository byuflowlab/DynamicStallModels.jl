using DelimitedFiles
using Plots
using FiniteDiff
using Statistics
using Polynomials

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end



polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/vertol 23010-1.58 experimental.csv", ',')
polar2 = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/xf-v23010-il-200000-n5.csv", ','; skipstart=12)

### Convert to radians
polar[:,1] = polar[:,1].*(pi/180)

tencl, tenidx = nearestto(polar[:,1], 10.0*pi/180)
linrng = [1, tenidx]

dcldalpha = gradient(polar[:,1], polar[:,2], polar[1:tenidx,1])
meandcldalpha = mean(dcldalpha)

# fun(x) = meandcldalpha*() #I'm not sure if this going to work... so 
# rng = (-pi, polar[end,1])
# alpha0 = find_zero(fun, rng, Bisection())

### I could just use y = mx + b.... for the linear range. and average that. => b = y - mx 
m = meandcldalpha
n = length(polar[:,1])
alpha0vec = zeros(n)
for i = 1:n
    local b
    b = polar[i,2] - m*polar[i,1]
    alpha0vec[i] = -b/m
end
meanalpha0 = mean(alpha0vec)

println("Average slope: ", meandcldalpha)
println("Average α_0: ", meanalpha0)
a = fit(polar[:,1], polar[:,2], 1)
b = fit(polar2[:,1].*(pi/180), polar2[:,2], 1)


plt = plot(polar[:,1], polar[:,2], legend=:bottomright, lab="Experimental")
plot!(polar2[:,1].*(pi/180), polar2[:,2], lab="Airfoil Tools")
plot!(a, lab="lsf exp")
plot!(b, lab="lsf af tools")
display(plt)


nothing


