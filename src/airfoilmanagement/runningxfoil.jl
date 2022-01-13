using Xfoil
using Plots
using Printf
using Polynomials

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

# extract geometry
x = Float64[]
y = Float64[]

f = open("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/coordinates/naca0012.dat", "r")

for line in eachline(f)
    entries = split(chomp(line))
    push!(x, parse(Float64, entries[1]))
    push!(y, parse(Float64, entries[2]))
end

close(f)

# set operating conditions
alpha = -15:1:30
re = 3.83e6

c_l, c_d, c_dp, c_m, converged = Xfoil.alpha_sweep(x, y, alpha, re, mach=0.383, iter=100, zeroinit=false, printdata=false)

aoa = alpha.*(pi/180) #Convert to radians

zeroalf, zeroidx = nearestto(c_l, 0.0)
tenalf, tenidx = nearestto(aoa, 10.0*pi/180)
linrng = zeroidx:tenidx

mxb = fit(aoa[linrng], c_l[linrng], 1)

dcldalpha = mxb.coeffs[2] #least squares fit

alpha0 = roots(mxb)[1] #least squares fit

inviscid(alf) = dcldalpha.*(alf.-alpha0)
cli = inviscid(aoa)


clplt = plot(alpha, c_l, legend=:bottomright)
plot!(alpha, cli)
display(clplt)

nothing