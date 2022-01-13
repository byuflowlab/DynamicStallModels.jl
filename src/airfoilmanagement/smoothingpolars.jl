using DelimitedFiles
using CCBlade
using FLOWMath
using Plots
using FiniteDiff
using Dierckx

include("/Users/adamcardoza/Box/research/FLOW/turbineopt/20kW_turbine/prework/src/airfoils/fixingairfoils.jl")

experimentaldata = "/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989staticdata_fig1.csv" #"/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/vertol 23010-1.58 experimental.csv"
polar = readdlm(experimentaldata, ',')

rads0 = polar[:,1].*(pi/180)
Cl0 = polar[:,2]
Cd0 = ones(length(polar[:,1])) #polar[:,3]
ar75 = 0.09602814324656855 *1 #The 75% chord/R_tip for the UAE, just for something to work against. 

rads, Cl, Cd = prefixdata(rads0, Cl0, Cd0; s=0.00)

# rd1, cl1, cd1 = rads, Cl, Cd

rads, Cl, Cd = viterna(rads, Cl, Cd, ar75)
# Cl[1] = 0.0 #Note: What am I doing here? -> It looks like I'm making sure the first entry is zero (duh), but I vaguely remember doing that because the Viterna method is supposed to? I forget why. 
# println(length(Cl))
Cl, Cd = smoothpolar(rads, Cl, Cd;nskip=1, peaks=true)
# println(length(Cl))

clfit = Spline1D(rads, Cl, k=2, s=0.25)
cl1 = clfit.(rads)

clmax, clmaxidx = findmax(Cl)
r1 = rads[clmaxidx]
r2 = -3*pi/180

slope = gradient(rads, Cl, rads)

slopefit = Spline1D(rads,slope,k=2,s=300)
slope1 = slopefit.(rads)
curvature = gradient(rads, slope1, rads)
jerk = gradient(rads, curvature, rads)

clplt = plot(rads, Cl, legend=:topleft, label="current")
plot!(rads, cl1, label="smoothed")
# scatter!(rads0, Cl0, label="static")
vline!([r1], label="α1")
vline!([r2], label="α2")
# xlims!((-10.0*pi/180, 20*pi/180))
display(clplt)

# slpplt = plot(rads,slope1, legend=false)
# vline!([r1])
# vline!([r2])

# crvplt = plot(rads, curvature, legend=false)
# vline!([r1])
# vline!([r2])
# ylims!((-200.0,200.0))

# jrkplt = plot(rads, jerk, legend=false)
# vline!([r1])
# vline!([r2])
# ylims!((-50.0,50.0))

# plt = plot(clplt, slpplt, crvplt, jrkplt, layout=(4,1))
# display(plt)

mat = hcat(rads, Cl, Cd)
# writedlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/extendedNACA0012_fromexperimental.dat", mat, ',')

nothing