using Plots, Statistics, DelimitedFiles, Roots, DynamicStallModels, DifferentialEquations


k = 0.2 #given
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

    amp = 2
    shift = 5

    alpha = amp*sin(omega*t) + shift
    return alpha*(pi/180)
end

function alphadot(t)
    # k = 0.1 #given
    # u = 50 #Assumed
    # c = 1.0 #Assumed

    # omega = k*2*u/c #implied

    amp = 2
    shift = 5

    alphadot = amp*omega*cos(omega*t)

    return alphadot*(pi/180)
end


aoa = -pi:0.01:pi
lift = 2*pi.*(aoa)
drag = zero(aoa)
polar = hcat(aoa, lift, drag)




A = [0.165, 0.335] #From the Hansen 2004 paper, for flat plate
b = [0.0455, 0.3000] # "" ""

Tp = 3.0 # "" "" 
Tf = 6.0 # "" "" 
T = [Tp, Tf]


m, n = size(polar)
newpolar = hcat(polar, zeros(m), zeros(m))
afs = [airfoil(newpolar; A, b, T)]

dsmodel = Riso(Functional(), length(afs), afs)

x0 = zeros(4)
x0[3] = 1.0

p = [U, Udot, alpha, alphadot, c]

tspan = (0.0, 80.0)


prob = ODEProblem(dsmodel, x0, tspan, p)
sol = solve(prob, dtmax=0.1)

Cl, Cd, t =  parsesolution(dsmodel, sol, p)



alphavec = alpha.(sol.t)

expdata = readdlm("./data/Hansen2004/figure8_flatplate/indicial.csv", ',')


clplt = plot(legend=:topleft, title="Cyclic Alpha", yaxis="Cl", xaxis="Alpha (deg)")
scatter!(expdata[:,1], expdata[:,2], lab="Hansen 2004")
plot!(alphavec.*(180/pi), Cl, lab="Riso")
plot!(alphavec.*(180/pi), linearlift.(alphavec), lab="Static")
display(clplt)



nothing