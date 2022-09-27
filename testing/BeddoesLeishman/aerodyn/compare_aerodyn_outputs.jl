using DynamicStallModels, OpenFASTsr, FLOWMath, Plots, Plots.PlotMeasures, DelimitedFiles

#=
Test the Beddoes-Leishman model as given in the documentation. I had to make some minor modifications as not all of the equations were given. 

Adam Cardoza 8/22/22
=#

of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)


# include("../../../src/BeddoesLeishman/BeddoesLeishmanAeroDyn.jl")
# include("../../../src/solve.jl")

# /Users/adamcardoza/.julia/dev/DynamicStallModels/data/
fullout = readdlm("../../../data/aerodynout_fordynamicstall.csv", ',')
# fullout = readdlm("../../../data/aerodyn_outputs.csv", ',')
names = fullout[1,:]
names[1] = "Time"
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in 1:length(names)) #11 or 12th node is what we're looking for. 

### Create airfoil
c = 1.0

du21_a17 = of.read_airfoilinput("../../../data/DU21_A17.dat")

dcndalpha = du21_a17.c_nalpha
alpha0 = du21_a17.alpha0*(pi/180)
polar = hcat(du21_a17.aoa.*(pi/180), du21_a17.cl, du21_a17.cd, du21_a17.cm)
clfit = Akima(polar[:,1], polar[:,2])
cdfit = Akima(polar[:,1], polar[:,3])
cmfit = Akima(polar[:,1], polar[:,4])

alphasep = [du21_a17.alpha2, du21_a17.alpha1].*(pi/180)
# alphasep = [-8, 14].*(pi/180) #Flow was seperating when I didn't expect it to. 

A = [du21_a17.a1, du21_a17.a2]
b = [du21_a17.b1, du21_a17.b2]
T = [du21_a17.t_p, du21_a17.t_f0, du21_a17.t_v0, du21_a17.t_vl]

S1 = du21_a17.s1
S2 = du21_a17.s2
S3 = du21_a17.s3
S4 = du21_a17.s4
S = [S1, S2, S3, S4]
xcp = 0.2

af = Airfoil(polar, clfit, cdfit, cmfit, dcndalpha, alpha0, alphasep, A, b, T, S, xcp)
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af




### Create model
zeta = -3 #Low-pass-filter frequency cutoff
A5 = du21_a17.a5
b5 = du21_a17.b5
Tsh = du21_a17.st_sh #Strouhal Frequency
eta = du21_a17.eta_e #Recovery factor
constants = [zeta, A5, b5, Tsh, eta]
dsmodel = BeddoesLeishman(Indicial(), 1, airfoils, 2, constants)


### Create simulation data
tvec = outs["Time"]
nt = length(tvec)

Uvec = [sqrt(outs["AB1N011Vx"][i]^2 + outs["AB1N011Vy"][i]^2) for i in 1:nt] #m/s
aoavec = outs["AB1N011Alpha"].*(pi/180)






### Solve
states, loads = solve_indicial(dsmodel, [c], tvec, Uvec, aoavec)

Cn = loads[:,1]



staticCn = clfit.(aoavec)

cnplt = plot(xaxis="time (s)", yaxis="Cn", right_margin=20mm, leg=:topleft)
plot!(tvec, Cn, lab="BL_AD")
plot!(tvec, staticCn, lab="static")
plot!(tvec, outs["AB1N011Cn"], linestyle=:dash, lab="OpenFAST")
plot!(twinx(), tvec, aoavec.*(180/pi), lab="AOA", leg=:bottomright, linecolor=:purple, linestyle=:dot, yaxis="AOA")
display(cnplt)

#=
Hey hey, that matches pretty stinkin well. That suggests that at least the model in the attached region is working properly. -> I'll need to get something that experiences a little more deflection and see how it compares. 
Adam Cardoza 8/25/22

=#

# nt, na = size(states)

# for i= 1:na
#     plt = plot(tvec, states[:,i], xaxis="Time (s)", leg=false, title="State $i" )
#     display(plt)
# end

nothing