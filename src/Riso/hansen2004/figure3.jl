using Plots, DelimitedFiles, FLOWMath, Roots

#=
Recreating figure 3 from Hansen's 2004 paper.

Seeking to make the seperation point function more robust. 

##### Problems
Todo. The Seperation point function has an odd dip between -20 and 0. That shouldn't be there. -> The dip is due to the value of dcldalpha. It doesn't behave nicely if the slope is too high. Question: But if the slopw is too low, then won't the inviscid lift under predict? 

Todo: The seperation point function looks like it states at 1 from about -10 degrees to about 7.5 degrees. Mine might be that way, but it looks like it isn't quite. 





###### Solved Problems. 

Todo. The seperation point function doesn't have any discontinuities in the paper. I have one at both the positive and negative angles of stall. -> I was artificially setting the function value to zero when the aoa was outside of afm and afp as suggested by the paper. The function natrually drives to zero.

#Todo. There is a jump in my fully seperated lift, whereas Hansen's doesn't have that. -> Due to the discontinuity in the seperation point function. 
=#

include("../riso.jl")

polar = readdlm("../../../experimentaldata/Hansen2004/figure3_separationfunction/static.csv", ',')

polar[:,1] = polar[:,1].*(pi/180)

liftfit = Akima(polar[:,1], polar[:,2])

alpha0 = find_zero(liftfit, 0.0)
dcldalpha = 2*pi*1.0
clmin, clmin_idx = findmin(polar[:,2])
clmax, clmax_idx = findmax(polar[:,2])
afm = polar[clmin_idx, 1] 
afp = polar[clmax_idx, 1]
linearlift(alpha) = dcldalpha*(alpha-alpha0)


alphas = [find_seperation_alpha(liftfit, dcldalpha, alpha0)...]
alfas = alphas.*(180/pi)

alfavec = collect(-50:0.5:50)
alphavec = alfavec.*(pi/180)
n = length(alphavec)

clvec = zeros(n)
clsvec = zeros(n)
clivec = zeros(n)
fvec = zeros(n)
for i = 1:n
    clvec[i] = Clfs(alphavec[i], liftfit, dcldalpha, alpha0) #Clfs(alphavec[i], liftfit, dcldalpha, alpha0, afm, afp)
    clsvec[i] = liftfit(alphavec[i])
    clivec[i] = linearlift(alphavec[i])
    fvec[i] = seperationpoint(alphavec[i], alphas[2], alphas[1], liftfit, dcldalpha, alpha0)
end








plt = plot(legend=:bottomright, xaxis="Angle of Attack (deg)", ylim=(-1.5, 1.5))
plot!(alfavec, clsvec, lab="Static Lift")
plot!(alfavec, clivec, lab="Inviscid Lift")
plot!(alfavec, clvec, lab="Fully Separated Lift")
plot!(alfavec, fvec, lab="Separation Point")
vline!(alfas, lab="Seperation AOA")
display(plt)


nothing