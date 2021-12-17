using Plots, DelimitedFiles, FLOWMath, Roots

"""
11/5/21 Adam Cardoza

Dr. Ning said that I should be able to completely recreate figure 3 from Hansen's 2004 paper. Exactly. So we'll do that, hopefully really quickly. 
"""

include("../Riso.jl")

polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure3_separationfunction/static.csv", ',')

polar[:,1] = polar[:,1].*(pi/180)

liftfit = Akima(polar[:,1], polar[:,2])

alpha0 = find_zero(liftfit, 0.0)
dcldalpha = 2*pi*1.02
linearlift(alpha) = dcldalpha*(alpha-alpha0)



alfavec = collect(-40:0.1:40)
alphavec = alfavec.*(pi/180)
n = length(alphavec)

clvec = zeros(n)
clsvec = zeros(n)
clivec = zeros(n)
fvec = zeros(n)
for i = 1:n
    clvec[i] = Clfs(alphavec[i], liftfit, dcldalpha, alpha0)
    clsvec[i] = liftfit(alphavec[i])
    clivec[i] = linearlift(alphavec[i])
    fvec[i] = fst(alphavec[i], liftfit, dcldalpha, alpha0)
end

plt = plot(legend=:topleft, xaxis="Angle of Attack (deg)", ylim=(-1.5, 1.5))
plot!(alfavec, clsvec, lab="Static Lift")
plot!(alfavec, clivec, lab="Inviscid Lift")
plot!(alfavec, clvec, lab="Fully Separated Lift")
plot!(alfavec, fvec, lab="Separation Point")
display(plt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure3_separationfunction/recreatedfigure3.png")

nothing