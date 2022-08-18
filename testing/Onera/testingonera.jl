using FLOWMath, DelimitedFiles, DifferentialEquations, Plots, Statistics


include("onera.jl")

function U(t)
    return 60
end


function alpha(t)
    k = 0.062
    v = 60
    c = 1.5

    amp = 4.85
    shift = 14.9

    omega = 2*k*v/c

    alfa =  amp*cos(omega*t) + shift
    return alfa*pi/180
end

function alphadot(t)
    k = 0.062
    v = 60
    c = 1.5

    amp = 4.85
    shift = 14.9

    omega = 2*k*v/c

    alfadot =  -amp*sin(omega*t)/omega
    return alfadot*pi/180
end

polar2 = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Larsen2007/radians.csv", ',') #readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Larsen2007/vertol_lowaoa_static.csv", ',') #Todo: Polar 1 is a actually a little off from the experimental data in Larsen's paper. 
# polar2[:,1] = polar2[:,1].*(pi/180)
# liftfit = Akima(polar[:,1], polar[:,2])
liftfit = Akima(polar2[:,1], polar2[:,2])
lfdvec = gradient.(Ref(polar2[:,1]), Ref(polar2[:,2]), polar2[:,1])
liftfitderiv = Akima(polar2[:,1], lfdvec)
dcldalpha = 6.653921403540334 #2*pi*1.05 #This has been checked. 

A = [0.3, 0.1]

omegahat = [0.125, 0.085] #TODO: I wonder if they rounded the values to "nice" values. 
omega = omegahat.*(2*60/1.5)
zeta = 0.7

p = [alpha, alphadot, liftfit, liftfitderiv, dcldalpha, A, omega, zeta]

tspan = (0.0, 5.0)
x0 = zeros(3)

prob = ODEProblem(states!,x0,tspan,p)
sol = solve(prob, dtmax=0.01)

statesplt = plot(sol,linewidth=2,xaxis="t",label=["x1" "x2" "x3" "x4"],layout=(4,1))
display(statesplt)


Cl1, t1, u1 = parsesolution(sol, p)
alphavec = alpha.(t1)
alfavec = alphavec.*(180/pi)

explift = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Larsen2007/dynamic_lift.csv", ',')

clplt = plot(leg=:bottomleft, title="Cyclic Lift")
plot!(alfavec, Cl1, lab="My ONERA")
scatter!(explift[:,1].*(180/pi), explift[:,2], lab="Paper Values")
display(clplt)

nothing
