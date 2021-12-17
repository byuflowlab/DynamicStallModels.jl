using Plots
using DelimitedFiles

function clinv(alpha; dcl=2*pi*1.27, alpha0=-0.048) #NOTE: I'm fudging the slope at zero lift value to match the experimental rather than calculating it
    return dcl*(alpha-alpha0)
end

function cl(alpha)
    polar = readdlm("/Users/adamcardoza/Box/research/FLOW/learning/exampledata/NACA4412.dat", '\t'; skipstart=3)
    clfit = Akima(polar[:,1], polar[:,2])
    #alphan, idxn = nearestto(polar[:,1], alpha)
    return clfit(alpha) #polar[idxn,2]
end

include("Oye.jl")

fs = 0.5 #fst #Assuming that we are starting at the static separation f, but we'll play
c = 1
Vrel = 15

aoavec = collect(-15:0.1:15) 
alphavec = aoavec.*(pi/180) #polar[:,1]
clvec = cl.(alphavec)
clinvvec = clinv.(alphavec)
clfsvec = clfs.(alphavec)
Clvec = Cl.(Ref(fs), alphavec, Ref(c), Ref(Vrel))

# clplt = plot(aoavec, clvec, lab="Static", leg=:topleft)
# plot!(aoavec, clinvvec, lab="Inviscid")
# plot!(aoavec, clfsvec, lab="Fully separated")
# plot!(aoavec, Clvec, lab="Dynamic")
# display(clplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/oye/steadyliftcurve.png")

# function alpha(t)
#     what = 0.001 #0.062
#     v = 60
#     c = 1.5
#     return 4.85*pi*cos(what*2*v*t/c)/180 + 9.0*pi/180
# end

function alpha(t) #Constant alpha
    return 8*sin(t)*pi/180
end

function U(t)
    return 60
end

tspan = (0.0, 300.0)
c = 1.5

Cld, tvec = Oye(tspan, alpha, U, c)
tvec2 = collect(0:15:300)
clst = cl.(alpha.(tvec2))

cldplt = plot(tvec, Cld, lab="dynamic") #, ylims=[0.8, 1.6]
scatter!(tvec2, clst, lab="static", markershape=:circle)
display(cldplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/oye/Cld_steadyaoa.png")