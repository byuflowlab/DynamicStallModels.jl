using DynamicStallModels, OpenFASTsr, FLOWMath, Plots, Plots.PlotMeasures, DelimitedFiles, LaTeXStrings, Revise

#=
Test the Beddoes-Leishman model as given in the documentation. I had to make some minor modifications as not all of the equations were given. 

Adam Cardoza 10/3/22
=#

of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)

include("parseaerodyn.jl")

# items = ["t", "alpha", "Cc_pot", "fprime_c", "fprimeprime_c", "Df", "Cn", "Cc", "PI"] 
items = ["t", "tau", "Cn_FS", "Cn_v", "Cn_inviscid", "Cn_alpha_q_nc", "Cn_q_circ", "Cn_alpha_q_circ", "f''", "Cn_alpha_nc", "Cn_q_nc", "T_alpha", "Kalpha_f", "Kprime_alpha", "M", "q_f_cur", "U", "c", "alpha", "Cn_viscous"]
file = "/Users/adamcardoza/.julia/dev/DynamicStallModels/data/aerodyn_intermediates.txt"

entries = readdlm(file, ',')

mat = parseaerodyn(entries, items)


# /Users/adamcardoza/.julia/dev/DynamicStallModels/data/
# fullout = readdlm("../../../data/gonzalezouts.csv", ',')
fullout = readdlm("../../../data/gonzalez2.csv", ',')
# fullout = readdlm("../../../data/aerodyn_outputs.csv", ',')
names = fullout[1,:]
names[1] = "Time"
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in 1:length(names)) #11 or 12th node is what we're looking for. 

### Create airfoil
c = 3.542
idx = 1

# du21_a17 = of.read_airfoilinput("../../../data/DU21_A17.dat")
du21_a17 = of.read_airfoilinput("../../../data/DU40_A17.dat")

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

# af = Airfoil(polar, clfit, cdfit, cmfit, dcndalpha, alpha0, alphasep, A, b, T, S, xcp)
af, constants = of.make_dsairfoil(du21_a17)
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af




### Create model
# zeta = -3 #Low-pass-filter frequency cutoff
# A5 = du21_a17.a5
# b5 = du21_a17.b5
# Tsh = du21_a17.st_sh #Strouhal Frequency
# eta = du21_a17.eta_e #Recovery factor
# constants = [zeta, A5, b5, Tsh, eta]
dsmodel = BeddoesLeishman(Indicial(), 1, airfoils, 3, constants)


### Create simulation data
tvec = outs["Time"]
nt = length(tvec)

if idx <10
    num = "00$idx"
elseif idx<100
    num ="0$idx"
else
    num = "$idx"
end

# Uvec = [sqrt(outs["AB1N011Vx"][i]^2 + outs["AB1N011Vy"][i]^2) for i in 1:nt] #m/s
Uvec = outs["AB1N"*num*"VRel"]
aoavec = outs["AB1N"*num*"Alpha"].*(pi/180)

Uvec[3:end] = mat[:,17]
aoavec[3:end] = mat[:,19]

errfun(x, xt) = (x-xt)/xt

# tvec = 0:0.1:5.0
# Uvec = Uvec[1].*ones(length(tvec))
# aoavec = aoavec[1].*ones(length(tvec))

### Solve
states, loads = solve_indicial(dsmodel, [c], tvec, Uvec, aoavec)

Cn = loads[:,1]



staticCn = clfit.(aoavec) #Todo: This will need to be rotated... methinks

cnplt = plot(xaxis="time (s)", yaxis="Cn", right_margin=20mm, leg=:topright)
plot!(tvec, Cn, lab="BL_AD")
plot!(tvec, staticCn, lab="static")
plot!(outs["Time"], outs["AB1N"*num*"Cn"], linestyle=:dash, lab="OpenFAST")
# plot!(twinx(), tvec, aoavec.*(180/pi), lab="AOA", leg=:bottomright, linecolor=:purple, linestyle=:dot, yaxis="AOA")
# display(cnplt)


staticCc = cdfit.(aoavec)
Cc = loads[:,2]
Cd0 = af.cd(af.alpha0)

ccplt = plot(xaxis="Time (s)", yaxis=L"C_c", right_margin=20mm, leg=:topright)
plot!(tvec, Cc, lab="BL_ADG") 
plot!(tvec, staticCc, lab="static")
plot!(outs["Time"], outs["AB1N011Ct"], linestyle=:dash, lab="OpenFAST")
# plot!(twinx(), tvec, aoavec.*(180/pi), lab="AOA", leg=:bottomright, linecolor=:purple, linestyle=:dot, yaxis="AOA")
# display(ccplt)

# Cterr = errfun.(Cc, outs["AB1N011Ct"])

#=
Most of the time it is hovering just below 1% error. (median is -0.72% error.) I'm going to guess this is due to the fact that I'm not using a tangential separation point function. It could also be due to the fact that I might need to rotate the slope of the lift curve. Although, the slope that I'm currently using from the AeroDyn data is the slope of the normal curve, so.... it should be solid. 
=#

# nt, na = size(states)

# for i= 1:na
#     plt = plot(tvec, states[:,i], xaxis="Time (s)", leg=false, title="State $i" )
#     display(plt)
# end



#=
12/9/22
The angles of attack that are getting passed to my solver and the ones getting passed to their solver are slightly different... I don't know if it's enough... I mean... it looks like its 0.018%.... so insignificant. So it's something else. 
=#

Cfsn, Cvn, Nnc_aq, Nc_q, Nc_aq, Nnc_a, Nnc_q, Talpha, M = extractintermediatestates(states, Uvec, c, af)

CNFSplt = plot(mat[:,1], mat[:,3], lab="OpenFAST", ylab="Cnfs")
plot!(tvec, Cfsn, lab="DSM")
# display(CNFSplt)

Cvnplt = plot(mat[:,1], mat[:,4], lab="OpenFAST", ylab="Cvn")
plot!(tvec, Cvn, lab="DSM")
# display(Cvnplt)

tauplt = plot(mat[:,1], mat[:,2], lab="OpenFAST", ylab="Tau")
plot!(tvec, states[:, 22], lab="DSM")
# display(tauplt)

#=
12/12/22
Okay number 1 problem, there is a significant discrepancy for CNFS. Then another discrepancy is that their tau value just keeps increasing (from a later starting point but it appears to be the same slope), but their Cvn value never changes, even though they have a tau value... Which... is weird. I should be getting a value here. 

Oddly... their CNFS plus their Cvn value does not equal the normal force value that they put out... which implies that they are returning a different Cn value... Their Cnfs value is fairly near their Cn value... It's making me think that they've flagged this node so that it always returns their non dynamic stall value. ... Ahh... but their Cnfs value matches their inviscid value, which when converted to their viscous value it matches.... so it is calculating there. 

While I was looking through, I found this UA_BlendSteady() function... that got me worried that perhaps that the point is getting the steady influence or something... but I don't think so, because it has a comment to that says to check if they should have turned off unsteady aerodynamics. 

I guess I just need to keep going down the rabbit hole til I find where things are different. 
=#

# Nnc_aqplt = plot(mat[:,1], mat[:,6], lab="OpenFAST", ylab="Nnc_aq")
# plot!(tvec, Nnc_aq, lab="DSM")
# display(Nnc_aqplt)

#Mine is constant, and theirs drops off.... 

# Ncirc_qplt = plot(mat[:,1], mat[:,7], lab="OpenFAST", ylab="Ncirc_q")
# plot!(tvec, Nc_q, lab="DSM")
# display(Ncirc_qplt)

#Same again, mine is constant, and theirs drops off. How odd....

# Nnc_aplt = plot(mat[:,1], mat[:,10], lab="OpenFAST", ylab="Nnc_a")
# plot!(tvec, Nnc_a, lab="DSM")
# display(Nnc_aplt)

#Again, mine is constant. That just doesn't make sense. Am I reading in the data incorrectly? 

# Talphaplt = plot(mat[:,1], mat[:,12], lab="OpenFAST", ylab="Talpha")
# plot!(tvec, Talpha, lab="DSM")
# display(Talphaplt)

# Both constant, -2.43% error, they are both very very small. 

Kaplt = plot(mat[:,1], mat[:,11], lab="OpenFAST", ylab="Ka")
plot!(tvec, states[:,4], lab="DSM")
display(Kaplt)

#=
I may have found the error. This is remaining zero, and OpenFAST does a little arc. 

#Changed to use the exact U and aoa that OpenFAST passes to the dynamic stall model. 

The difference in how I do the initialization makes me out of phase and have a larger initial amplitude. 
=# 

Kpaplt = plot(mat[:,1], mat[:,12], lab="OpenFAST", ylab="Kpa")
plot!(tvec, states[:,10], lab="DSM")
display(Kpaplt) #Todo: There is still a constant difference here. 



qplt = plot(mat[:,1], mat[:,16], lab="OpenFAST", ylab="q", leg=:bottomright)
plot!(tvec, states[:,3], lab="DSM")
display(qplt)

#=
Interesting... our low pass constants are giving different values. #Todo. Currently comparing low pass constants. But I need to be looking at why our pitching values are different. 

There appears to be slight differences in the input velocity. I'm going to guess there are slight differences in the input angle as well... They're probably rounding for output. -> Their filter cutoff value is actually 0.5.... so I don't know where it gets changed to that... but that is really odd. #Todo: Find where the filter cutoff value is changed. 

Alright, Now I'm passing in the exact U and aoa values that OpenFAST is passing to the dynamic stall model and now q converges to the values that they have. 


=#

nothing