using DelimitedFiles
using Plots
using FLOWMath

include("BeddoesLeishman.jl")


function U(t)
    return 60
end

function Udot(t)
    return 0
end


function alpha(t)
    what = 0.062
    v = 60
    c = 1.5
    amp = 4.85
    shift = 7.33 #14.92 #9.0
    conv = pi/180 
    return amp*conv*cos(what*2*v*t/c)/180 + shift*conv
end

# function alphadot(t)
#     what = 0.062
#     v = 60
#     c = 1.5
#     amp = 4.85
#     return -what*2*v*t*amp*pi*sin(what*2*v*t/c)/(180*c)
# end

# function CL0(t)
#     what = 0.062
#     v = 60
#     c = 1.5
#     amp = 4.85
#     shift = 7.33 #14.92 #9.0

#     alpha = amp*pi*cos(what*2*v*t/c)/180 + shift*pi/180
#     alpha0 = -0.019
#     dcldalpha = 2*pi
#     return dcldalpha*(alpha-alpha0)
# end

# function CL0dot(t)
#     what = 0.062
#     v = 60
#     c = 1.5
#     dcldalpha = 2*pi*1.05
#     amp = 4.85
    
#     return -dcldalpha*what*2*v*t*amp*pi*sin(what*2*v*t/c)/(180*c)
# end


#Enviroment
V = 60

#Geometry
c = 1.5
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/learning/exampledata/xf-v23010-il-200000-n5.csv", ','; skipstart=12)
# liftfit = Akima(polar[:,1].*(pi/180), polar[:,2])
polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/extendedVertol 23010-1.58_fromexperimental.dat", ',')
liftfit = Akima(polar[:,1], polar[:,2])
dcldalpha = 2*pi
alpha0 =0 # -0.019



#Finding separation angles of attack
aoa = polar[:,1]
maxcl, maxidx = findmax(polar[:,2])
alphav = aoa[maxidx] #14.75*pi/180 #What the report said was pretty much the same as the static stall angle. 
# plot(aoa, polar[:,2], leg=false)
# vline!(aoa[[maxidx]])


v = 60

#Constants
A = [0.3, 0.7, 1.0, 1.0]
b = [0.0455, 0.3] #Larsen   #[0.14, 0.53] #Leishman89 
# Tp = 1/0.4125
# Tf = 1/0.0875
S = [0.05, 0.05]
tspan = [0.0, 1.5]
alpha1 = alphav

Cld, tvec = BeddoesLeishman(tspan, alpha, U, c, A, b, S, alpha0, alpha1, dcldalpha;Tp=1.7, Tf=3.0, Tv=6.0, nt=100, f0=NaN, a=343.0)



alphavec = (alpha.(tvec)).*(180/pi)
Cls = liftfit.(alpha.(tvec))
aoavec = collect(-0.05:0.05:12*pi/180)
clo =  dcldalpha.*(aoavec.-alpha0)  # CL0.(tvec)
expdata = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Larsen2007/larsenexperimental_970a.csv", ',')

# clplt = plot(tvec, Cld, lab="dynamic", leg=:topright)

# plot!(tvec, Cls, lab="static") #TODO: note that the static data only goes to around 15 degrees and so as the aoa oscillates, it is oscillating out of the region of provided data. 
# xlabel!("time (s)")
# ylabel!("Coefficient of Lift")
# display(clplt)
# # savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishman/Cld_stepaoa.png")

cycleplot = plot(alphavec, Cld, lab="dynamic", leg=false)
plot!(alphavec, Cls, lab="static")
plot!(aoavec.*(180/pi), clo, lab="inviscid")
scatter!(expdata[:,1], expdata[:,2], lab="experimental")
xlabel!("angle (degrees)")
ylabel!("Coefficient of Lift")
display(cycleplot)

# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishman/clvaoa71321a.png")

nothing