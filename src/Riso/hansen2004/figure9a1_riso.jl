using Plots, Statistics, DelimitedFiles, Roots

include("../Riso.jl")

k = 0.1 #given
u = 1.0 #Assumed
c = 1.0 #Assumed

omega = k*2*u/c #implied

function U(t)
    return u
end

function Udot(t)
    return 0.0
end

function V(t)
    return 0.0
end

function alpha(t)
    # k = 0.1 #given
    # u = 50 #Assumed
    # c = 1.0 #Assumed

    # omega = k*2*u/c #implied

    amp = 4
    shift = 3

    alpha = amp*sin(omega*t) + shift
    return alpha*(pi/180)
end

function alphadot(t)
    # k = 0.1 #given
    # u = 50 #Assumed
    # c = 1.0 #Assumed

    # omega = k*2*u/c #implied

    amp = 4
    shift = 3

    alphadot = amp*omega*cos(omega*t)

    return alphadot*(pi/180)
end

polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/static.csv", ',')

plr = deepcopy(polar)
plr[:,1] = plr[:,1].*(pi/180)
liftfit = Akima(plr[:,1], plr[:,2])

dcldalpha = 2*pi*1.05
alpha0 = find_zero(liftfit, 0.0)

yintercept = -alpha0*dcldalpha

linearlift(alpha) = dcldalpha*alpha + yintercept
linearcl = linearlift.(plr[:,1])

staticplt = plot(legend=:topleft, title="Static Lift Plot") #, xlim=(-0.07,0.25))
scatter!(plr[:,1], plr[:,2], lab="Static")
plot!(plr[:,1], linearcl, lab="Linear")
# display(staticplt)


A = [0.294, 0.331] #From the Hansen 2004 paper, for Riso airfoil
b = [0.0664, 0.3266] # "" ""
Tp = 1/0.4125 #From the Larsen 2007 paper for the Vertol airfoil
Tf = 1/0.0875 # "" "" 
# A = [ 0.3, 0.7] #Leishman 1990
# b = [0.14, 0.53] # "" "" 
# Tp = 3.0 # "" "" 
# Tf = 6.0 # "" "" 



x0 = zeros(4)
p = [c, A, b, Tp, Tf, dcldalpha, liftfit, U, Udot, V, alpha, alphadot, alpha0]
tspan = (0.0, 80.0)

prob = ODEProblem(states!, x0, tspan, p)
sol = solve(prob)

Cld1, Cdd1, Cmd1, u1, t1 = parsesolution(sol, p)
alphavec = alpha.(sol.t)

# tspan2 = (0.0, 3.0)

# x02 = u1[end,:]
# prob2 = ODEProblem(states!, x02, tspan2, p)
# sol2 = solve(prob2)

# Cld2, u2, f2 = parsesolution(sol2, p)
# alphavec = alpha.(sol2.t)

# statesplt = plot(sol2,linewidth=2,xaxis="t",label=["x1" "x2" "x3" "x4"],layout=(4,1))
# display(statesplt)


expdata = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/circ3pm4.csv", ',')


clplt = plot(legend=:topleft, title="Cyclic Alpha", yaxis="Cl", xaxis="Alpha (deg)")
scatter!(expdata[:,1], expdata[:,2], lab="Paper Values")
plot!(alphavec.*(180/pi), Cld1, lab="My Values")
plot!(alphavec.*(180/pi), liftfit.(alphavec), lab="Static")
display(clplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/riso/nonlinearsolve/attachedlift.png")


nothing