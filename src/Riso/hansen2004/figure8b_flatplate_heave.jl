using Plots, Statistics, DelimitedFiles, Roots

include("../Riso.jl")

k = 0.2 #given
u = 1.0 #Assumed
udot = 0.0
c = 1.0 #Assumed
theta = 5.0 #given
thetadot = 0.0

omega = k*2*u/c #implied

function Y(t)
    return 0.2*c*sin(omega*t)
end

function V(t)
    return 0.2*c*omega*cos(omega*t)
end

function Vdot(t)
    return -0.2*c*omega*omega*sin(omega*t)
end

function U(t)
    return sqrt(u^2 + V(t)^2)
end

function Udot(t)
    return (2*u*udot - 2*V(t)*Vdot(t))/(2*sqrt(u^2 + V(t)^2))
end

function alpha(t)
    phi = atan(V(t), u)


    alpha = phi + theta*(pi/180)
    return alpha
end

function alphadot(t)
    phidot = (u*Vdot(t) - V(t)*udot)/((u^2)*(((V(t)/u)^2)+1))

    alphadot = phidot + thetadot*(pi/180)

    return alphadot
end

# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure8_flatplate/static.csv", ',')
polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure8_flatplate/flatplatepolar.csv", ',')

plr = deepcopy(polar)
plr[:,1] = plr[:,1].*(pi/180)
liftfit = Akima(plr[:,1], plr[:,2])

dcldalpha = 2*pi
alpha0 = 0.0 #find_zero(liftfit, 0.0)

yintercept = -alpha0*dcldalpha

linearlift(alpha) = dcldalpha*alpha + yintercept
linearcl = linearlift.(plr[:,1])

staticplt = plot(legend=:topleft, title="Static Lift Plot") #, xlim=(-0.07,0.25))
scatter!(plr[:,1], plr[:,2], lab="Static")
plot!(plr[:,1], linearcl, lab="Linear")
# display(staticplt)


A = [0.165, 0.335] #From the Hansen 2004 paper, for flat plate
b = [0.0455, 0.3000] # "" ""
# A = [ 0.3, 0.7] #Leishman 1990
# b = [0.14, 0.53] # "" "" 
Tp = 3.0 # "" "" 
Tf = 6.0 # "" "" 



x0 = zeros(4)
p = [c, A, b, Tp, Tf, dcldalpha, liftfit, U, Udot, V, alpha, alphadot, alpha0]
tspan = (0.0, 80.0)

prob = ODEProblem(states!, x0, tspan, p)
sol = solve(prob, dtmax=0.1)

Cld1, u1, f1 = parsesolution(sol, p)
yc = Y.(sol.t)./c


expdata = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure8_flatplate/heave.csv", ',')


clplt = plot(legend=:topleft, title="Cyclic Alpha", yaxis="Cl", xaxis="y/c", ylims=(0.0, 1.1))
scatter!(expdata[:,1], expdata[:,2], lab="Paper Values")
plot!(yc, Cld1, lab="My Values")
plot!(yc, linearlift.(alpha.(sol.t)), lab="static")
display(clplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/riso/flatplate/recreatefig8b.png")



nothing