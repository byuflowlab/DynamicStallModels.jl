using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations

dsm = DynamicStallModels


path = dirname(@__FILE__)
cd(path)

file = "../../polars/Extended_NACA_0030.csv"

polar = readdlm(file , ',')

c = 0.55
M = 0.11
v = 343
Vrel = M*v

dsmodel = Riso(Functional(), [0.294, 0.331] , [0.0664, 0.3266], [0.010933, 0.043732]) #add a way for varying Tp and Tf values
af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.RSP())

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af

tspan = (0.0, 2.0)

function Uvector(t)
    return 0.11*343
end

function alphavec(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.0
    k = 0.051

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function alphavecdot(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.0
    k = 0.051

    v = M*a
    omega = k*2*v/c

    alf = amp*omega*cos(omega*t)
    return alf*(pi/180)
end

parameters = [Uvector, 0.0, alphavec, alphavecdot]
x_initial = [0.0, 0.0, 0.0, 0.0]

prob = ODEProblem(airfoils, x_initial, tspan, parameters)
sol = DifferentialEquations.solve(prob, reltol=1e-8)

answer, extra, stuff, more = parsesolution(dsmodel, airfoils, sol, parameters)

plot(answer[1,:], extra[1,:])
plot!(answer[1,:], stuff[1,:])
plot!(answer[1,:] , answer[2,:])
