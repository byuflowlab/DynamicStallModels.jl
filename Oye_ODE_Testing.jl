using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, OpenFASTsr, LaTeXStrings, DifferentialEquations

dsm = DynamicStallModels
of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)

#original c = 0.55, M = 0.07

c = 0.55

M = 0.11
a = 343.0
Vrel = M*a #60

af = airfoil(polar_0015; A=[8], sfun=dsm.LSP()) 



airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af

dsmodel = Oye(Functional(), 1, airfoils, 1, 2)

function Uvector(t)
    return Vrel
end

function alpha(t)
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

p = [[4], c, Uvector, alpha]

tspan = (0.0 , 2.0)

x0 = [0.7]

prob = ODEProblem(dsmodel, x0, tspan, p)

sol = DifferentialEquations.solve(prob)

Lift, aoa = parsesolution_Oye(sol, p, af)

plot(aoa, Lift)