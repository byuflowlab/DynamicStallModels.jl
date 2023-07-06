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

dsmodel = Riso(Indicial(), [0.294, 0.331] , [0.0664, 0.3266], [1.5, 6]) 
af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.RSP())

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af

tvec = range(0.0, 2.0, 1000)

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

avec = alphavec.(tvec)
uvec = Vrel.*ones(length(tvec))
avecdot = alphavecdot.(tvec)

states, loads = solve_indicial(airfoils, tvec, uvec, avec; alphadotvec = avecdot)

plot(avec, loads[:,1])