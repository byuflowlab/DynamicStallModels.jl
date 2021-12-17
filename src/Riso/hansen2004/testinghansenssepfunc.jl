using Plots, FLOWMath, DelimitedFiles

"""
11/1/21 Adam Cardoza
Trying to figure out how to deal with Cl diverging to infinity as f->1. Also trying to establish the behavior of f that I expect (high in the linear region and low past that.) I think... I should look up what it is supposed to be. 

"""

include("../Riso.jl")

# function fs(alpha, liftfit, dcldalpha, alpha0)
#     return ((2*sqrt(liftfit(alpha)/(dcldalpha*(alpha-alpha0))))-1)^2
# end

# function clfs(alpha, liftfit, dcldalpha, alpha0)
#     f = fs(alpha, liftfit, dcldalpha, alpha0)
#     return (liftfit(alpha)- (dcldalpha*(alpha-alpha0)*f))/(1-f)
# end

polar1 = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/extendedVertol 23010-1.58.dat", ',')
polar2 = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Larsen2007/vertol_lowaoa_static.csv", ',')
polar2[:,1] = polar2[:,1].*(pi/180)
# liftfit = Akima(polar[:,1], polar[:,2])
liftfit = Akima(polar2[:,1], polar2[:,2])

dcldalpha = 2*pi*1.05
alpha0 = 0.014 # -0.019 #This is definitely negative. 

clifun(alpha) = dcldalpha*(alpha-alpha0)

alfavec = collect(2.0:0.1:12.8)
alphavec = alfavec.*(pi/180)
n = length(alphavec)

clvec = zeros(n)
clsvec = zeros(n)
clivec = zeros(n)
for i = 1:n
    clvec[i] = Clfs(alphavec[i], liftfit, dcldalpha, alpha0)
    clsvec[i] = liftfit(alphavec[i])
    clivec[i] = clifun(alphavec[i])
end

clplt = plot(;leg=:topleft, xlim=(0,15), xaxis="Angle of Attack (deg)", yaxis="Lift Coefficient")
plot!(alfavec, clvec, lab="Fully Separated")
plot!(alfavec, clsvec, lab="Static")
plot!(alfavec, clivec, lab="Inviscid")
scatter!(polar2[:,1].*(180/pi), polar2[:,2], lab="exp from paper")
# scatter!(polar1[:,1].*(180/pi), polar1[:,2], lab="static_exp")
display(clplt)

fvec = zeros(length(alphavec))
for i = 1:length(alphavec)
    fvec[i] = fst(alphavec[i], liftfit, dcldalpha, alpha0)
end

fplt = plot(;leg=:topleft, title="f")
plot!(alfavec, fvec)
# display(fplt)
nothing