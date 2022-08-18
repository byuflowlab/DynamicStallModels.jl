using DelimitedFiles
using Plots
using FLOWMath
using Polynomials

include("../Larsen.jl")


function U(t)
    M = 0.379
    a = 343.0
    return M*a
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
    c = 0.1
    M = 0.379
    a = 343.0
    shift = 10.3
    amp = 8.1
    k = 0.075

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function alphadot(t)
    c = 0.1
    M = 0.379
    a = 343.0
    shift = 10.3
    amp = 8.1
    k = 0.075

    v = M*a
    omega = k*2*v/c
    
    alfd = amp*omega*cos(omega*t)
    return alfd*(pi/180)
end

function alphaddot(t)
    c = 0.1
    M = 0.379
    a = 343.0
    shift = 10.3
    amp = 8.1
    k = 0.075

    v = M*a
    omega = k*2*v/c
    
    alfdd = -amp*(omega^2)*sin(omega*t)
    return alfdd*(pi/180)
end

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

### Polar Analysis
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/xf-n0012-il-1000000.csv", ','; skipstart=12)
polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989staticdata_fig1.csv", ',')


polar[:,1] = polar[:,1].*(pi/180) #Convert to radians

liftfit = Akima(polar[:,1], polar[:,2])

zeroalf, zeroidx = nearestto(polar[:,2], 0.0) #Nearest to zero lift
tenalf, tenidx = nearestto(polar[:,1], 10.0*pi/180) #Ten is typically before stall #*pi/180
linrng = zeroidx:tenidx

mxb = fit(polar[linrng,1], polar[linrng,2], 1)

dcldalpha = mxb.coeffs[2] #least squares fit 

alpha0 = roots(mxb)[1] #least squares fit   # NACA 0012 is a symmetric airfoil... so shouldn't the zero lift be at zero? 

# polarplt = plot(polar[:,1], polar[:,2], leg=false)
# plot!(mxb)
# display(polarplt)

maxcl, maxclidx = findmax(polar[:,2]) #critical lift coefficient (given by Larsen 2007), the beginning of LE separation. 
# println("Static Max Cl: ", maxcl)
alphav = polar[maxclidx, 1] 

maxcl = dcldalpha*(alphav-alpha0) #Projecting the inviscid lift to the static stall point
maxcl = 1.2

# function CL0(t)
#     c = 0.1
#     M = 0.379
#     a = 343.0
#     shift = 10.3
#     amp = 8.1
#     k = 0.075

#     v = M*a
#     omega = k*2*v/c

#     alpha = shift + amp*sin(omega*t)
#     # alpha0 = -0.019
#     # dcldalpha = 2*pi
#     return dcldalpha*(alpha-alpha0)
# end

# function CL0dot(t)
#     c = 0.1
#     M = 0.379
#     a = 343.0
#     shift = 10.3
#     amp = 8.1
#     k = 0.075

#     v = M*a
#     omega = k*2*v/c
    
#     alfdd = -amp*(omega^2)*sin(omega*t)
    
#     return dcldalpha*alfdd*(pi/180)
# end


### Constants
M = 0.379
a = 343.0
v = M*a #Putting the speed here to calculate w. 

### Geometry
c = 0.1

#Constants #Todo: I'm not sure if these are the correct constants. I think A will be the same. and probably b. and the time constants tend to stay the same from iteration to iteration. However... they might be different because of the speed difference. 
A = [0.165, 0.335].*.7
# b = [0.0455, 0.3]
# Tp = 1/0.4125
# Tf = 1/0.0875

#Finding separation angles of attack
# aoa = polar[:,1]
# maxcl, maxidx = findmax(polar[:,2])
# alphav = aoa[maxidx] #14.75*pi/180 #What the report said was pretty much the same as the static stall angle. 
# plot(aoa, polar[:,2], leg=false)
# vline!(aoa[[maxidx]])


w = ones(4)
w[1] = 0.0455*2*v/c
w[2] = 0.3*2*v/c
w[3] = 0.1*2*v/c
w[4] = 0.075*2*v/c



#Initialize 
u0 = zeros(4)
u0[1] = 0.0 #Fudging the initial conditions to have a derivative of zero initially. Which produced much better agreement. 
u0[2] = 0.06
u0[3] = 1.0
u0[4] = 0.1
p = [CL0, CL0dot, alpha, alphadot, alphav, w]
tspan = (0.0, 0.05)

prob = ODEProblem(states!, u0, tspan, p)
sol = solve(prob, dtmax=0.0001)

statesplt = plot(sol,linewidth=2,xaxis="t",title=["x1" "x2" "θd" "Clv"],layout=(4,1), leg=false)
# display(statesplt)


Cld = parsesolution(sol, p)
alphavec = alpha.(sol.t)
Cls = liftfit.(alphavec)
# aoavec = collect(-0.05:0.05:12*pi/180)
# clo =  dcldalpha.*(aoavec.-alpha0)  # CL0.(sol.t)
exppolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989fig8Cn.csv", ',')

alfavec = alphavec.*(180/pi)
mat = hcat(alfavec, Cld)
# writedlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/outputs/leishman1989/naca0012_8.1_10.3s_Cdn_Larsen.csv", mat, ',')

# clplt = plot(sol.t, Cld, lab="dynamic", leg=:topright)
# plot!(sol.t, Cls, lab="static") #TODO: note that the static data only goes to around 15 degrees and so as the aoa oscillates, it is oscillating out of the region of provided data. 
# xlabel!("time (s)")
# ylabel!("Coefficient of Lift")
# display(clplt)
# # savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/larsen/Cld_stepaoa.png")

cycleplot = plot(alpha.(sol.t).*(180/pi), Cld, lab="Larsen", leg=:topleft)
# plot!(alpha.(sol.t).*(180/pi), Cls, lab="static")
# plot!(aoavec.*(180/pi), clo, lab="inviscid")
scatter!(exppolar[:,1], exppolar[:,2], lab="experimental")
xlabel!("angle (degrees)")
ylabel!("Coefficient of Lift")
display(cycleplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/larsen/larsensclvaoa71321a.png")