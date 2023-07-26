using DynamicStallModels, OpenFASTsr, FLOWMath, Plots, Plots.PlotMeasures, DelimitedFiles, LaTeXStrings, Revise, Statistics

#=
Test the Beddoes-Leishman model as given in the documentation. I had to make some minor modifications as not all of the equations were given. 

Adam Cardoza 10/3/22
=#

#Todo: alpha (generic function with 5 methods). I'm not sure that there should be 5 methods. 

of = OpenFASTsr
DSM = DynamicStallModels

path = dirname(@__FILE__)
cd(path)

include("parseaerodyn.jl")

# items = ["t", "alpha", "Cc_pot", "fprime_c", "fprimeprime_c", "Df", "Cn", "Cc", "PI"] 
items = ["i", "t", "U", "alpha", "a_s", "M", "Gonzalez_factor", "c", "C_nalpha", "A1", "b1", "A2", "b2", "eta_e", "tau_v", "Cn_FS", "Cn_v", "Cn_alpha_q_nc", "Cn_q_circ", "Cn_alpha_q_circ", "fprimeprime", "Cn_alpha_nc", "Cn_q_nc", "T_alpha", "Kalpha_f", "Kprime_alpha", "q_f_cur", "alpha_filt_cur", "k_alpha", "alpha_e", "Kq_f", "Kprime_q", "Df", "fprime", "alpha_f", "Cc_pot", "fprime_c", "fprimeprime_c", "C_nalpha_circ", "Cn", "Cc", "Cm"] #, "Cn_visc", "Cc_visc"]
file = "/Users/adamcardoza/.julia/dev/DynamicStallModels/data/ADG_intermediates.txt"

entries = readdlm(file, ',')

mat = parseaerodyn(entries, items, 19)


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

du21_a17 = of.read_airfoilinput("../../../data/airfoils/DU40_A17.dat")

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
af = of.make_dsairfoil(du21_a17, c)
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af




### Create model
# zeta = -3 #Low-pass-filter frequency cutoff
# A5 = du21_a17.a5
# b5 = du21_a17.b5
# Tsh = du21_a17.st_sh #Strouhal Frequency
# eta = du21_a17.eta_e #Recovery factor
# constants = [zeta, A5, b5, Tsh, eta]
# dsmodel = BeddoesLeishman(Indicial(), 1, airfoils, 3)


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
aoavec = (outs["AB1N"*num*"Alpha"].+0.01334).*(pi/180)

Uvec[3:end] = mat[:,3, idx]
aoavec[3:end] = mat[:,4, idx]

errfun(x, xt) = (x-xt)/xt

# tvec = 0:0.1:5.0
# Uvec = Uvec[1].*ones(length(tvec))
# aoavec = aoavec[1].*ones(length(tvec))

### Solve
states, loads = solve_indicial(airfoils, tvec, Uvec, aoavec)

Cn = zero(tvec)
Cc = zero(tvec)

for i in eachindex(tvec)
    Cn[i], Cc[i] = DSM.rotate_load(loads[i,1], loads[i,2], states[i, 1])
end 




cnplt = plot(xaxis="time (s)", yaxis="Cn", right_margin=20mm, leg=:bottomright)
plot!(tvec, Cn, lab="BL_AD")
# plot!(outs["Time"], outs["AB1N"*num*"Cn"], linestyle=:dash, lab="OpenFAST")
plot!(mat[:,2, idx], mat[:,40, idx], lab="OpenFAST intermediate", linestyle=:dash)
# plot!(twinx(), tvec, aoavec.*(180/pi), lab="AOA", leg=:bottomright, linecolor=:purple, linestyle=:dot, yaxis="AOA")
# display(cnplt)
# # savefig(cnplt, "OpenFASTcomparison_cn_modifiedNREL5MW.png")

# cnerror = errfun.(Cn[3:end], mat[:,40]).*100
#-> Now there is less than half a percent error. Which.... I'm quite happy with. Like.... it could/should be better... but... it seems good enough. Their output seems heavily rounded. Like there is no oscillation. Which means that my solver should account for more fatique? or just more noise... one of the two. -> When comparing to the intermediate step (the easily accessed output is rounded down), the max percent relative error is 0.183%. Which... is pretty good methinks. I don't know what is causing the difference? 


Cd0 = af.cd(af.alpha0)

ccplt = plot(xaxis="Time (s)", yaxis=L"C_c", right_margin=20mm, leg=:bottomright)
plot!(tvec, Cc, lab="DSM") 
# plot!(outs["Time"], outs["AB1N001Ct"], linestyle=:dash, lab="OpenFAST")
plot!(mat[:,2, idx], mat[:,41, idx], lab="OpenFAST intermediate", linestyle=:dash)
# plot!(twinx(), tvec, aoavec.*(180/pi), lab="AOA", leg=:bottomright, linecolor=:purple, linestyle=:dot, yaxis="AOA")
# display(ccplt)
# # savefig(ccplt, "OpenFASTcomparison_cc_modifiedNREL5MW.png")


# Cterr = errfun.(Cc, outs["AB1N001Ct"]).*100
# ccerror = errfun.(Cc[3:end], mat[:,41]).*100


#=
Max error is 0.8%, which is pretty good, but not as low as I'd like. I'd like to say that it's negligible because the total difference is pretty small. Like the difference we'd see in the total loading should be pretty small, right? Well... if they're both scaled by the same factor, we should still be 0.8% off... so it still has a higher error than I'd like. But... this is the rounded error, so it might be less. 

Well.. that brought it down to 0.71%. Still not quite as close as I'd like. I really wonder what the difference is. -> I was rotating the final time by alpha instead of aoa.. which decreased error to 0.65%. 
=#

Cm = loads[:,3]

cmplt = plot(xaxis="Time (s)", yaxis=L"C_m", right_margin=20mm, leg=:bottomright)
plot!(tvec[3:end], Cm[3:end], lab="DSM") 
plot!(mat[:,2,idx], mat[:,42,idx], lab="OpenFAST")
# display(cmplt)
# # savefig(cmplt, "OpenFASTcomparison_cm_modifiedNREL5MW.png")





# nt, na = size(states)

# for i= 1:na
#     plt = plot(tvec, states[:,i], xaxis="Time (s)", leg=false, title="State $i" )
#     display(plt)
# end



#=
12/9/22
The angles of attack that are getting passed to my solver and the ones getting passed to their solver are slightly different... I don't know if it's enough... I mean... it looks like its 0.018%.... so insignificant. So it's something else.  -> 12/22/22 Well... actually...  it was making a difference. 
=#

##############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# Cfsn, Cvn, Nnc_aq, Nc_q, Nc_aq, Nnc_a, Nnc_q, Talpha, M, k_alpha, TI, alphae, k_q, Cpot, dcndalpha_circ = DSM.extractintermediatestates(states, Uvec, c, af; a=335.0)

# TI_of = mat[:,18]./mat[:,22] #c/a #The same. 

# CNFSplt = plot(mat[:,1], mat[:,3], lab="OpenFAST", ylab="Cnfs")
# plot!(tvec, Cfsn, lab="DSM")
# # display(CNFSplt) 
# # Todo. Constant and off. -> I wonder if the Cvn is actually getting applied. -> The fp was off, so it was throwing this off. #Todo: There is this large hit at the beginning, but other than that, it looks like it is spot on. 

# Cvnplt = plot(mat[:,1], mat[:,4], lab="OpenFAST", ylab="Cvn")
# plot!(tvec, Cvn, lab="DSM")
# # display(Cvnplt)
# #Todo: Constant for most of it, then mine picks up but theirs doesn't. The overall variance isn't large.... so I don't think it's a problem. 

# tauplt = plot(mat[:,1], mat[:,2], lab="OpenFAST", ylab="Tau")
# plot!(tvec, states[:, 22], lab="DSM")
# # display(tauplt)
# #Todo: Theirs constantly increases, mine increases and drops repeatedly. They appear to be the same slope. 

# #=
# 12/12/22
# Okay number 1 problem, there is a significant discrepancy for CNFS. Then another discrepancy is that their tau value just keeps increasing (from a later starting point but it appears to be the same slope), but their Cvn value never changes, even though they have a tau value... Which... is weird. I should be getting a value here. 

# Oddly... their CNFS plus their Cvn value does not equal the normal force value that they put out... which implies that they are returning a different Cn value... Their Cnfs value is fairly near their Cn value... It's making me think that they've flagged this node so that it always returns their non dynamic stall value. ... Ahh... but their Cnfs value matches their inviscid value, which when converted to their viscous value it matches.... so it is calculating there. 

# While I was looking through, I found this UA_BlendSteady() function... that got me worried that perhaps that the point is getting the steady influence or something... but I don't think so, because it has a comment to that says to check if they should have turned off unsteady aerodynamics. 

# I guess I just need to keep going down the rabbit hole til I find where things are different. 
# =#

# Nc_aq_plt = plot(mat[:,1], mat[:,8], lab="OpenFAST", ylab=L"N^c_{aq}")
# plot!(tvec, Nc_aq, lab="DSM")
# # display(Nc_aq_plt)  #Good. 

# Nnc_aqplt = plot(mat[:,1], mat[:,6], lab="OpenFAST", ylab="Nnc_aq")
# plot!(tvec, Nnc_aq, lab="DSM")
# # display(Nnc_aqplt)
# #Todo. This is off. It almost looks like it is phase shifted. ... The values are tiny... so I don't know if it'll even have that much of an effect on the overall outcome. I'm just worried that the difference will make a difference later done the line. -> Fixed this components. 

# Cnoncirc_a_plt = plot(mat[:,1], mat[:,29], lab="OpenFAST", ylab=L"N^{nc}_\alpha")
# plot!(tvec, Nnc_a, lab="DSM") 
# # display(Cnoncirc_a_plt) #Good

# Cnc_q_plt = plot(mat[:,1], mat[:,30], lab="OpenFAST", ylab=L"N^{nc}_q")
# plot!(tvec, Nnc_q, lab="DSM") 
# # display(Cnc_q_plt) #Todo. Oscillates the wrong way. Fixed. I was missing a negative. 

# Kqplt = plot(mat[:,1], mat[:,31], lab="OpenFAST", ylab=L"K_q")
# plot!(tvec, states[:,5], lab="DSM")
# # display(Kqplt) #Todo. Off from each other. I bet there's something wrong with the filtering here as well. ... Well... I shifted from using the filtered values to using the unfiltered values, like they did... but that made things worse. -> Added q_f as a state, now it matches. 

# Kpq_plt = plot(mat[:,1], mat[:,32], lab="OpenFAST", ylab=L"K'_q")
# plot!(tvec, states[:,11], lab="DSM")
# # display(Kpq_plt) #TODO: Follows decently well, there is some odd discontinuity after the first pass but then quickly converges. 





# Ncirc_qplt = plot(mat[:,1], mat[:,7], lab="OpenFAST", ylab="Ncirc_q")
# plot!(tvec, Nc_q, lab="DSM")
# # display(Ncirc_qplt) #Todo. This appears phase shifted. -> I hadn't filtered pitch rate to using 

# alphaeplt = plot(mat[:,1], mat[:,28], lab="OpenFAST", ylab=L"\alpha_e")
# plot!(tvec, alphae, lab="DSM") 
# # display(alphaeplt)



# Nnc_aplt = plot(mat[:,1], mat[:,10], lab="OpenFAST", ylab="Nnc_a")
# plot!(tvec, Nnc_a, lab="DSM")
# # display(Nnc_aplt)

# #Again, mine is constant. That just doesn't make sense. Am I reading in the data incorrectly? 

# k_alpha_plt = plot(mat[:,1], mat[:,21], lab="OpenFAST", ylab=L"k_\alpha")
# plot!(tvec, k_alpha, lab="DSM")
# # display(k_alpha_plt) #It appears there is some fluctuation in the value, and it isn't constant. dcndalpha, M, A1, A2, b1, b2 all appear to be the same. An error of -0.0571%... I don't expect there to be any error, because that can compound... but I think I'm just going to have to accept it for now. 

# Talphaplt = plot(mat[:,1], mat[:,12], lab="OpenFAST", ylab="Talpha")
# plot!(tvec, Talpha, lab="DSM")
# # display(Talphaplt)
# # Todo. Both constant, -2.43% error, they are both very very small. -> Now it's the same error of -0.057%. 

# Kaplt = plot(mat[:,1], mat[:,13], lab="OpenFAST", ylab="Ka")
# plot!(tvec, states[:,4], lab="DSM")
# # display(Kaplt) #Todo. There's a difference here. -> They were updating Ka using the filtered q values after. 

# #=
# I may have found the error. This is remaining zero, and OpenFAST does a little arc. 

# #Changed to use the exact U and aoa that OpenFAST passes to the dynamic stall model. 

# The difference in how I do the initialization makes me out of phase and have a larger initial amplitude. 
# =# 

# Kpaplt = plot(mat[:,1], mat[:,14], lab="OpenFAST", ylab="Kpa", leg=:bottomright)
# plot!(tvec, states[:,10], lab="DSM")
# # display(Kpaplt) #Todo. There is still a constant difference here. -> I was plotting the wrong thing. There still is a difference, but I was plotting the wrong thing. This is probably off because Talpha is off, because it's a function of Ka and Talpha. Ka is on now, so Talpha must be off. -> There was a slight error in my extraction function, now it should be good. 



# qfplt = plot(mat[:,1], mat[:,16], lab="OpenFAST", ylab=L"q_f", leg=:bottomright)
# plot!(tvec, states[:,29], lab="DSM")
# # display(qfplt)

# alphaplt = plot(mat[:,1], mat[:,20], lab="OpenFAST", ylab=L"\alpha_{filt}", leg=:bottomright)
# plot!(tvec, states[:,1], lab="DSM")
# # display(alphaplt)

# #=
# Interesting... our low pass constants are giving different values. #Todo. Currently comparing low pass constants. But I need to be looking at why our pitching values are different. 

# There appears to be slight differences in the input velocity. I'm going to guess there are slight differences in the input angle as well... They're probably rounding for output. -> Their filter cutoff value is actually 0.5.... so I don't know where it gets changed to that... but that is really odd. #Todo. Find where the filter cutoff value is changed. -> I copied what they did in their files.... so we good. I guess. 

# Alright, Now I'm passing in the exact U and aoa values that OpenFAST is passing to the dynamic stall model and now q converges to the values that they have. 
# -> Turns out that I wasn't. The aoa values that are passed out of OpenFAST are rounded degrees, which apparently makes a difference. Now I've got aoa matching. Which made q match. 

# =#


# fppplt = plot(mat[:,1], mat[:,9], lab="OpenFAST", ylab=L"f''", leg=:bottomright, markershape=:x)
# plot!(tvec[3:end], states[3:end,19], lab="DSM")
# # display(fppplt) #Todo. Offset -> I fix fp.

# fpplt = plot(mat[:,1], mat[:,34], lab="OpenFAST", ylab=L"f'", leg=:bottomright, markershape=:x)
# plot!(tvec, states[:,18], lab="DSM")
# # display(fpplt) #Todo. Offset -> Used there separation point function. Like I almost copied and pasted it. They calculate fst at every iteration... which is probably faster. TODO: There is a blip at the beginning. OpenFAST seems to just cut it off. I think I'm just going to leave it for now and make the system just report the static cl and cd (and therefore the static cn and cc). 

# fperr = errfun.(states[3:end, 18], mat[:,34]).*100
# ### Todo: I think that fp is the source of the error in my Cn.... I wonder if there is a difference in the alpha_f

# alphafplt = plot(mat[:,1], mat[:,35], lab="OpenFAST", ylab=L"\alpha_f", leg=:bottomright)
# plot!(tvec, states[:,2], lab="DSM")
# # display(alphafplt) #Okay, this is good, which means that it's my separation point function. 

# alphaferr = errfun.(states[3:end,2], mat[:,35]).*100
# #The relative error of alpha_f is really low, the average is 4.06e-6%... so miniscule. 

# #=
# 12/29/22 Debugging Cc

# I'm going to guess that there is a problem with fppc and now that I'm looking at it, I see that they use a Ccpot... I have Cpot... which if I'm following my notation scheme, then it would be the same value, but there's always possibility for error there. 


# Eta-e checks out, they're both set to 1.0.

# The gonzalez factor also checks out. (both are 0.2)

# =#

# Cpotplt = plot(mat[:,1], mat[:,37], lab="OpenFAST", ylab=L"C_{pot}", leg=:topright)
# plot!(tvec, Cpot, lab="DSM")
# # display(Cpotplt) #Todo. Close, but phase shifted. -> They use the unfiltered alpha here.... which is problematic in my POV. It means that I have to add another state. -> Boom, it was the fact that they use the unfiltered aoa. 

# cpoterr = errfun.(Cpot[3:end], mat[:,37]).*100

# fppcplt = plot(mat[:,1], mat[:,38], lab="OpenFAST", ylab=L"f''_c", leg=:topright, markershape=:x)
# plot!(tvec, states[:,21], lab="DSM")
# # display(fppcplt) #Todo. There is a minor shift. -> But it might be significant. Fixed it. They were using the circulation dcndalpha. 

# fppcerr = errfun.(states[3:end,21], mat[:,38]).*100
# ### Todo: The error here is 0.47% on average... which isn't huge, but explains the error. W


# fpcplt = plot(mat[:,1], mat[:,40], lab="OpenFAST", ylab=L"f'_c", leg=:topright, markershape=:x)
# plot!(tvec, states[:,20], lab="DSM")
# # display(fpcplt)

# fpcerr = errfun.(states[3:end, 20], mat[:,40]).*100
# ### Todo: The error is the same as fppc, so it suggests that the delay function doesn't significantly increase the error. 

# dcndalpha_circplt = plot(mat[:,1], mat[:,41], lab="OpenFAST", ylab=L"\frac{dC_n}{d\alpha}_c", leg=:topright, markershape=:x)
# plot!(tvec, dcndalpha_circ, lab="DSM")
# # display(dcndalpha_circplt)

# dcnerr = errfun.(dcndalpha_circ[3:end], mat[:,41]).*100
# ### The mean error here is -5.93e-15%. Which is incredibly small. So I can safely assume that this is not wrong... at all. 




loadsplt = plot(xaxis="time (s)", yaxis="Cn", right_margin=20mm, leg=:outerright)
plot!(tvec, Cn, lab=L"DSM $C_n$")
plot!(tvec, Cc, lab=L"C_t") 
plot!(tvec[3:end], Cm[3:end], lab=L"C_m") 
plot!(mat[:,2,idx], mat[:,40,idx], lab=L"OpenFAST $C_n$", linestyle=:dash)
plot!(mat[:,2,idx], mat[:,41,idx], lab=L"C_t", linestyle=:dash)
plot!(mat[:,2,idx], mat[:,42,idx], lab=L"C_m", linestyle=:dash)
display(loadsplt)

nothing