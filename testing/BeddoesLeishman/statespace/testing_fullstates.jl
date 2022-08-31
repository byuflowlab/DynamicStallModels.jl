using Plots, DifferentialEquations, FLOWMath, DelimitedFiles, Polynomials, Statistics, LaTeXStrings

#=
Test the full 13 (14 if you include tau) state Beddoes-Leishman model. Compare against the NACA 0012 data. 

Adam Cardoza 8/24/22s

=#

include("../../../src/BeddoesLeishman/BeddoesLeishman.jl")

### Constants
mach = 0.3
a = 343.3
v = mach*a #Putting the speed here to calculate w. 

A = [0.3, 0.7, 1.5, -0.5]
b = [0.14, 0.53, 0.25, 0.1, 0.5] #Todo: Need new values here. 
S = [2.75, 1.4] 
K = [0.0175, -0.12, 0.04]
m = 2.0

T_p = 1.7
T_f = 3.0
T_v = 6.0
T_vl = 7.5
a = 343.3

dCndalpha = 0.113*180/pi
alpha1 = 8.0*(pi/180) #14 #Changing this as a large effect on the maximum lift and the minimum lift. -> TODO: This should probably be defined by where f=0.7, not by the stall on the static plot. 
alpha0 = 0.17*(pi/180)
Cn1 = 1 #1.31 #Todo: I find it oddly strange that the lift curve isn't an input. I'm positive that it's an input. 
Cmo = -0.0037

### Initialize 
u0 = zeros(13)
u0[1] = 0.00
u0[2] = 0.00
u0[3] = 0.00
u0[4] = 5e-6
u0[5] = 0.00
u0[6] = 0.00
u0[7] = 0.00
u0[8] = 0.00
u0[9] = 0.00
u0[10] = 1.00 #Assuming the vortex starts attached
u0[11] = 0.00 #Assuming the vortex starts attached
u0[12] = 0.00
u0[13] = 1.00

tspan = (0.0, 0.5)


kay = 0.10 #If kay=0.01, the model responds too quickly. -> Decreasing the reduced frequency does increase the drop, but the amount of time that the airfoil is stalled is low. 

dtminn = 0.0001
alg = Tsit5()

fvec = [4, 5.33, 6, 8, 10]
kcvec = [0.051, 0.067, 0.076, 0.102, 0.127]

ceevec = zeros(length(fvec))

for i = 1:length(fvec)
    ceevec[i] = kcvec[i]*mach*a/(fvec[i]*pi)
end

cee = round(mean(ceevec); sigdigits=4)

Ufun1, Udotfun1, alphafun1, alphadotfun1, alphaddotfun1 = prepenvironment(; c=cee, M=mach, a=343.3, amp=10.0, shift=5.0, k=kay)

p1 = [Ufun1, Udotfun1, alphafun1, alphadotfun1, cee, dCndalpha, alpha1, Cn1, Cmo, A, b, S, K, m, T_p, T_f, T_vl, T_v, a]

prob1 = ODEProblem(fullstates!, u0, tspan, p1)
sol1 = solve(prob1, alg; dtmin=dtminn, force_dtmin=true) #;dtmax=0.0001 solve(prob)

t1, u1, du1, Cnt1, Cdt1, Cmt1 = fullstates_parsesolution(sol1, p1; eta=0.97)

alphavec1 = alphafun1.(t1)
alfavec1 = alphavec1.*(180/pi)

println("Run Information: ")
println("Chord Length: ", cee)
println("")



polar = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/experimentaldata/Leishman1989/leishman1989staticdata_fig1.csv", ',')
polar[:,1] = polar[:,1].*(pi/180) #Convert to radians

attachedlift = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/experimentaldata/leishman1989statespace/fig8/naca0012/leishman1989 fig 8 lift naca0012.csv", ',')
mal = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/experimentaldata/leishman1989statespace/fig8/naca0012/leishman1989 fig 8 naca 0012 lift model.csv", ',')

attacheddrag = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/experimentaldata/leishman1989statespace/fig8/naca0012/Fig8_drag_exp.csv", ',')
mad = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/experimentaldata/leishman1989statespace/fig8/naca0012/fig8_drag_model.csv", ',')

attachedmnt = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/experimentaldata/leishman1989statespace/fig8/naca0012/fig8_moment_exp.csv", ',')
mam = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/experimentaldata/leishman1989statespace/fig8/naca0012/fig8_moment_model.csv", ',')


liftfit = Akima(polar[:,1], polar[:,2])

Csn1 = liftfit.(alphavec1)



cv, cvdot, ccn, ccndot, tstartvec, tstopvec, tvlvec, Hcv = analyzesolution(sol1, p1; eta=0.97)

ccndot2 = gradient(t1, ccn, t1)
cvdot2 = gradient(t1, cv, t1)



cyclepltattached = plot(legend=:topleft, xaxis="Angle of Attack (degrees)", yaxis="Normal Coefficient", title="NACA 0012 - Mach = $mach")
scatter!(attachedlift[:,1], attachedlift[:,2], lab="Experimental")
scatter!(mal[:,1], mal[:,2], lab="Model - paper")
plot!(alfavec1, Cnt1, lab="1st Model")
plot!(alfavec1, Csn1, lab="Static")
vline!([alpha1].*(180/pi), lab="Stall Angle")
display(cyclepltattached)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/comparison/compareblmodels_naca0012_10.0a_10.0s09162021_attachedflow.png")
#=
It appears that my model isn't responding in time so it doesn't drop. It looks like it might begin to respond, it just never does. Maybe the frequency is much too high? 
 
-> I wonder if I just solve it with a simple RK4... what would happen? 
=#

dragplt = plot(;legend=:topleft, xaxis="Angle of Attack (degrees)", yaxis="Drag Coefficient")
plot!(alfavec1, Cdt1, label="Beddoes Leishman")
scatter!(attacheddrag[:,1], attacheddrag[:,2], lab="Experimental")
scatter!(mad[:,1], mad[:,2], lab="Model - paper")
display(dragplt)

momentplt = plot(;legend=:topleft, xaxis="Angle of Attack (degrees)", yaxis="Moment Coefficient")
plot!(alfavec1, Cmt1, label="Beddoes Leishman")
scatter!(attachedmnt[:,1], attachedmnt[:,2], lab="Experimental")
scatter!(mam[:,1], mam[:,2], lab="Model - paper")
display(momentplt)

mat1 = hcat(u1, alfavec1)
plt1 = plot(t1, mat1[:,1], title="x1")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt2 = plot(t1, mat1[:,2], title="x2")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt3 = plot(t1, mat1[:,3], title="x3")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt4 = plot(t1, mat1[:,4], title="x4")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt5 = plot(t1, mat1[:,5], title="x5")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt6 = plot(t1, mat1[:,6], title="x6")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt7 = plot(t1, mat1[:,7], title="x7")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt8 = plot(t1, mat1[:,8], title="x8")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt9 = plot(t1, mat1[:,9], title="C\'n")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt10 = plot(t1, mat1[:,10], title="f\'\'")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt11 = plot(t1, mat1[:,11], title="tau")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt12 = plot(t1, mat1[:,12], title="Cvn")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt13 = plot(t1, mat1[:,13], title="fqs")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
plt14 = plot(t1, mat1[:,14], title="alpha")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
stateplt = plot(plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8, plt9, plt10, plt11, plt12, plt13, plt14, layout=14, leg=false, size=(1000,1000))
display(stateplt)

mat2 = hcat(du1, alfavec1)
pltd1 = plot(t1, mat2[:,1], title="dx1")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd2 = plot(t1, mat2[:,2], title="dx2")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd3 = plot(t1, mat2[:,3], title="dx3")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd4 = plot(t1, mat2[:,4], title="dx4")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd5 = plot(t1, mat2[:,5], title="dx5")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd6 = plot(t1, mat2[:,6], title="dx6")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd7 = plot(t1, mat2[:,7], title="dx7")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd8 = plot(t1, mat2[:,8], title="dx8")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd9 = plot(t1, mat2[:,9], title="dC\'n")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd10 = plot(t1, mat2[:,10], title="df\'\'")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd11 = plot(t1, mat2[:,11], title="dtau")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd12 = plot(t1, mat2[:,12], title="dCvn")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd13 = plot(t1, mat2[:,13], title="dfqs")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
pltd14 = plot(t1, mat2[:,14], title="alpha")
vline!(tstopvec, la=0.5)
vline!(tstartvec, la=0.5)
derivstateplt = plot(pltd1, pltd2, pltd3, pltd4, pltd5, pltd6, pltd7, pltd8, pltd9, pltd10, pltd11, pltd12, pltd13, pltd14, layout=14, leg=false, size=(1000,1000))
display(derivstateplt)

cvdotplt = plot(legend=:outerright,xaxis="time (s)", yaxis="Cvdot")
plot!(t1, cvdot2, lab="Numerical",markershape=:circle)
plot!(t1, cvdot, lab="Analytical")
vline!(tstartvec, la=0.5, lab="start")
vline!(tstopvec, la=0.5, lab="stop")
display(cvdotplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/comparison/cvderivatives.png")

cvplt = plot(;legend=:outerright, xaxis="time (s)", yaxis=L"C_v")
plot!(t1, cv, lab=L"C_v")
vline!(tstopvec, la=0.5, lab="stop")
vline!(tstartvec, la=0.5, lab="start")
# vline!(tvlvec, la=0.5, lab=L"2T_{vl}")
display(cvplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/comparison/cvplot.png")

ccndotplt = plot(legend=:outerright,xaxis="time (s)", yaxis="Ccndot")
plot!(t1, ccndot2, lab="Numerical",markershape=:circle)
plot!(t1, ccndot, lab="Analytical")
vline!(tstartvec, la=0.5, lab=false)
vline!(tstopvec, la=0.5, lab=false)
display(ccndotplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/comparison/ccnderivatives.png")

Ccnplt = plot(legend=:false, xaxis="Time (s)", yaxis=L"C^c_N")
plot!(t1, ccn)
vline!(tstopvec, la=0.5, lab=false)
vline!(tstartvec, la=0.5, lab=false)
display(Ccnplt)

hplt = plot(leg=:outerright)
plot!(t1, Hcv, lab=L"H_{Cv}")
vline!(tstopvec, la=0.5, lab="stop")
vline!(tstartvec, la=0.5, lab="start")
display(hplt)
#=
This plot currently tells me that the airfoil appears to always be in stall... but the Cv plot says otherwise. 

I thought that at first that the f value might be driving Cv down, but... 1) f never reaches 0.0 (Cv does), and 2) ohhh...
It might be driving to zero when f goes to 1. Does it go to 1? It does indeed, right where Cv drives to zero. Which makes sense by its equation. 

Interestingly, if I remember correctly, this bad boy should only trigger under a certain condition... and I'm wondering if Heavi should be 1 or zero under attached conditions. I'm thinking that it should be 0, but that's just my thought. Uhh... Yeah... I think I have my Heaviside function backwards. No no.... I have it right. If x[11] (tau) is less than 2*T_vl, than the argument of the Heavi will be positive, and will return 1. If it is greater than, then it will return 0. Which is what I expect. 


=#

nothing