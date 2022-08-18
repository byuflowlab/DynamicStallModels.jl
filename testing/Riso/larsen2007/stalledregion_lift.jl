using DelimitedFiles, Plots, Statistics

include("../../Riso.jl")

expdata = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Larsen2007/stalled_region_experimental.csv", ',')

middle = (maximum(expdata[:,1])+minimum(expdata[:,1]))/2
delalpha = maximum(expdata[:,1])-middle


function U(t)
    return 60
end

function Udot(t)
    return 0
end

function V(t)
    return 0.0
end

function alpha(t)
    k = 0.062
    v = 60
    c = 1.5

    amp = 4.81
    shift = 14.9

    omega = 2*k*v/c
    omega = k
    # omega = 2*k*v/(10*c)

    alfa =  amp*sin(omega*t) + shift
    return alfa*pi/180
end

function alphadot(t)
    k = 0.062
    v = 60
    c = 1.5

    amp = 4.81
    shift = 14.9

    omega = 2*k*v/c
    omega = k # I get a shape much more like Larsen's data when I use the reduced frequency instead of the frequency. 
    # Not divided by 10. That.... gives you a wonky shape. 
    # omega = 2*k*v/(10*c)

    alfadot =  amp*cos(omega*t)/omega
    return alfadot*pi/180
end

# function alpha(t) #step
#     if t<150
#         return 4.0*pi/180
#     else
#         return 15.0*pi/180
#     end
# end

# function alphadot(t)
#     if t<149.5
#         return 0.0
#     elseif 149.5<t<150.5
#         return 22.0*pi/180
#     else
#         return 0.0
#     end
# end

function alpha34(t) #Assuming that the aoa at the three-quarters point is the same as the aoa for now, but it really shouldn't be. 
    return alpha(t)
end

#Enviroment
v = 60

#Geometry
c = 1.5
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/learning/exampledata/xf-v23010-il-200000-n5.csv", ','; skipstart=12)
# liftfit = Akima(polar[:,1].*(pi/180), polar[:,2])
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/extendedVertol 23010-1.58.dat", ',')
polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Larsen2007/fullpolar.csv", ',')  

polar = hcat(polar, zeros(length(polar[:,1])), zeros(length(polar[:,1])))

# liftfit = Akima(polar[:,1], polar[:,2])
liftfit = Akima(polar[:,1], polar[:,2])
dcldalpha = 2*pi*1.05
alpha0 = 0.014 # -0.019

#Constants
A = [0.3, 0.7] #[0.165, 0.335] # I certified these values. 
b = [0.14, 0.53] #[0.0455, 0.3] #[0.455, 0.3]
Tp = 1.7 #6.0 #1/0.4125
Tf = 11 #1/0.0875

#Finding separation angles of attack
# alphas = findsepalpha(liftfit, dcldalpha, alpha0)
# aoa = polar[:,1].*(pi/180)
# cli = dcldalpha.*(polar[:,1].-alpha0) 
# cliplt = plot(aoa, polar[:,2], leg=false)
# plot!(aoa, cli)
# display(cliplt)


#Initialize 
x0 = zeros(4)
x0[1] = 0.03
x0[2] = 0.1
x0[3] = 1.6
x0[4] = 0.34
p = [c, A, b, Tp, Tf, dcldalpha, liftfit, U, Udot, V, alpha, alphadot, alpha0]
tspan = (0.0, 1800.0)

prob = ODEProblem(states!,x0,tspan,p)
sol = solve(prob, dtmax=0.05)

statesplt = plot(sol,linewidth=2,xaxis="t",label=["x1" "x2" "x3" "x4"],layout=(4,1))
display(statesplt)


Cld, Cdd, Cmd, u, t = parsesolution(sol, p, polar)
alphavec = alpha.(sol.t)
Cls = liftfit.(alphavec)
alfavec = alphavec.*(180/pi)

clplt = plot(sol.t, Cld, lab="dynamic", leg=:bottomright)
plot!(sol.t, Cls, lab="static") #TODO: note that the static data only goes to around 15 degrees and so as the aoa oscillates, it is oscillating out of the region of provided data. 
xlabel!("time (s)")
ylabel!("Coefficient of Lift")
# display(clplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/riso/Cld_stepaoa.png")

# fplt = plot(sol.t, f, leg=false)
# xlabel!("time (s)")
# ylabel!("Separation factor")
# display(fplt)

staticcl = liftfit.(alphavec)

cycleplt = plot(legend=:topleft)
plot!(alfavec, Cld, lab="Riso")
scatter!(expdata[:,1], expdata[:,2], lab="Experimental")
plot!(alfavec, staticcl, lab="Static")
display(cycleplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/riso/larsen/attachedregion_lift.png")


nothing