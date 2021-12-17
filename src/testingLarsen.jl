using DelimitedFiles
using Plots
using FLOWMath

include("Larsen.jl")


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
    return amp*pi*cos(what*2*v*t/c)/180 + shift*pi/180
end

function alphadot(t)
    what = 0.062
    v = 60
    c = 1.5
    amp = 4.85
    return -what*2*v*amp*pi*sin(what*2*v*t/c)/(180*c)
end

function CL0(t)
    what = 0.062
    v = 60
    c = 1.5
    amp = 4.85
    shift = 7.33 #14.92 #9.0

    alpha = amp*pi*cos(what*2*v*t/c)/180 + shift*pi/180
    alpha0 = -0.019
    dcldalpha = 2*pi
    return dcldalpha*(alpha-alpha0)
end

function CL0dot(t)
    what = 0.062
    v = 60
    c = 1.5
    dcldalpha = 2*pi*1.05
    amp = 4.85
    
    return -dcldalpha*what*2*v*amp*pi*sin(what*2*v*t/c)/(180*c)
end


#Enviroment
V = 60

#Geometry
c = 1.5
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/learning/exampledata/xf-v23010-il-200000-n5.csv", ','; skipstart=12)
# liftfit = Akima(polar[:,1].*(pi/180), polar[:,2])
polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/extendedVertol 23010-1.58_fromexperimental.dat", ',')
liftfit = Akima(polar[:,1], polar[:,2])
dcldalpha = 2*pi
alpha0 = -0.019


#Constants
A = [0.165, 0.335]
b = [0.0455, 0.3]
Tp = 1/0.4125
Tf = 1/0.0875

#Finding separation angles of attack
aoa = polar[:,1]
maxcl, maxidx = findmax(polar[:,2])
alphav = aoa[maxidx] #14.75*pi/180 #What the report said was pretty much the same as the static stall angle. 
# plot(aoa, polar[:,2], leg=false)
# vline!(aoa[[maxidx]])

# c = 1
v = 60
w = ones(4)
w[1] = 0.0455*2*v/c
w[2] = 0.3*2*v/c
w[3] = 0.1*2*v/c
w[4] = 0.075*2*v/c



#Initialize 
u0 = zeros(4)
u0[1] = -0.015 #Fudging the initial conditions to have a derivative of zero initially. Which produced much better agreement. 
u0[2] = -0.002
u0[3] = 1.05
p = [CL0, CL0dot, alpha, alphadot, alphav, w]
tspan = (0.0, 1.5)

prob = ODEProblem(states!, u0, tspan, p)
sol = solve(prob, dtmax=0.0001)

# statesplt = plot(sol,linewidth=2,xaxis="t",label=["x1" "x2" "x3" "x4"],layout=(4,1))
# display(statesplt)


Cld = parsesolution(sol, p)
alphavec = alpha.(sol.t)
Cls = liftfit.(alphavec)
aoavec = collect(-0.05:0.05:12*pi/180)
clo =  dcldalpha.*(aoavec.-alpha0)  # CL0.(sol.t)
expdata = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Larsen2007/larsenexperimental_970a.csv", ',')

# clplt = plot(sol.t, Cld, lab="dynamic", leg=:topright)
# plot!(sol.t, Cls, lab="static") #TODO: note that the static data only goes to around 15 degrees and so as the aoa oscillates, it is oscillating out of the region of provided data. 
# xlabel!("time (s)")
# ylabel!("Coefficient of Lift")
# display(clplt)
# # savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/larsen/Cld_stepaoa.png")

cycleplot = plot(alpha.(sol.t).*(180/pi), Cld, lab="Larsen", leg=:topleft)
plot!(alpha.(sol.t).*(180/pi), Cls, lab="static")
plot!(aoavec.*(180/pi), clo, lab="inviscid")
scatter!(expdata[:,1], expdata[:,2], lab="experimental")
xlabel!("angle (degrees)")
ylabel!("Coefficient of Lift")
display(cycleplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/larsen/larsensclvaoa71321a.png")