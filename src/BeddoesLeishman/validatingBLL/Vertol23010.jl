using DelimitedFiles
using Plots
using FLOWMath

include("../BeddoesLeishmanLarsen.jl")


function U(t)
    return 60
end

function Udot(t)
    return 0
end

function V(t)
    return 0
end

function Vdot(t)
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

function alphaddot(t)
    what = 0.062
    v = 60
    c = 1.5
    amp = 4.85
    return -what^2*4*v^2*amp*pi*cos(what*2*v*t/c)/(180*c^2)
end



#Geometry
c = 1.5
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/learning/exampledata/xf-v23010-il-200000-n5.csv", ','; skipstart=12)
# liftfit = Akima(polar[:,1].*(pi/180), polar[:,2])
polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/extendedVertol 23010-1.58_fromexperimental.dat", ',')
liftfit = Akima(polar[:,1], polar[:,2])
dcldalpha = 6.369467381507933*1.05 #lsf
#6.7219 #average
# 6.369467381507933 #lsf  
#2*pi #guess 
#Makes a slight change in the slope of the inviscid solution, which consequently affects the upper region of the solution. -> I got this from doing a linear regression of the slop of the static data, and then extending that linearly to zero. 
alpha0 = 0.009 #Playing     
#0.01146683 #lsf
#0.01579 #Average  
#0.009 #Playing  
aoa = polar[:,1]
maxcl, maxidx = findmax(polar[:,2])
maxcl = 1.6 #critical lift coefficient (given by Larsen 2007), the beginning of LE separation.
alphav = .244 # 
#aoa[maxidx] # = 0.2262426801152533 
#14.75*pi/180 
#The angle that corresponds to f = 0.7, which closely corresponds to the static stall angle of attack (from Leishman 1989, just after equation 17)

### Constants
v = 60 #Putting the speed here to calculate w. 
w = ones(7)
w[1] = 0.125*2*v/c
w[2] = 0.375*2*v/c
w[3] = 0.275*2*v/c
w[4] = 0.075*2*v/c
w[5] = 2.5*2*v/c
w[6] = 2.5*2*v/c
w[7] = 0.4*2*v/c

A = [0.3, 0.7, 1.0, 1.0]
b = [0.0455, 0.3]
S = [0.05, 0.05]
a = 343.0

### Initialize 
u0 = zeros(8)
u0[1] = 0.0  
u0[2] = 0.0 
u0[3] = 0.00
u0[4] = 0.0 
u0[5] = 0.0
u0[6] = 0.75
u0[7] = 0.00 #Assuming the vortex starts attached
u0[8] = 0.00

states!, p, u0 = BeddoesLeishman(U, Udot, V, Vdot, alpha, alphadot, alphaddot, c, dcldalpha, alpha0, alphav, maxcl, A, S, w; u0=u0)
tspan = (0.0, 3.0)

### Solve the system
prob = ODEProblem(states!, u0, tspan, p)
sol = solve(prob)

# statesplt = plot(sol,linewidth=2,xaxis="t",title=["c1" "c2" "c3" "c4" "cl0d" "fd" "tau" "Clv"], leg=false,layout=(2,4))
# display(statesplt)


Cld, Cdd = parsesolution(sol, p)
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

cycleplot = plot(alpha.(sol.t).*(180/pi), Cld, lab="BLL", leg=:topleft)
# plot!(alpha.(sol.t).*(180/pi), Cls, lab="static")
# plot!(aoavec.*(180/pi), clo, lab="inviscid")
scatter!(expdata[:,1], expdata[:,2], lab="experimental")
xlabel!("angle (degrees)")
ylabel!("Coefficient of Lift")
ylims!((0.0,1.5))
display(cycleplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanlarsen/larsen2007fig8/initialbllclvaoa080421_matched.png")
