using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations

dsm = DynamicStallModels
# of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)
include("./NACA_0015_DSM.jl")

c = 0.55

M = 0.11
a = 343.0
Vrel = M*a #60

# polar = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/learning/exploring/exampledata/NACA4412.dat", '\t'; skipstart=3) 
#polar = readdlm("./data/polars/naca0012.txt", skipstart=3)
#polar[:,1] = polar[:,1].*pi/180
# polar_0015 = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/learning/exploring/exampledata/NACA0015.dat", '\t'; skipstart=3)
polar = polar_0015

dsmodel = Oye(Functional(), 1, 2, 4.0)

#du21_a17 = of.read_airfoilinput("../../data/airfoils/DU40_A17.dat") 
#af = of.make_dsairfoil(du21_a17, c) 

af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.LSP())

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af



#Note: alphasep is much higher for Faber's implemenation of the dsmodel. -> It might need more tuning... but it's something. 


tspan = (0, 2.0) #0:0.001:0.05

function Uvector(t)
    return 0.11*343.0
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

parameters = [Uvector, 0.0, alpha, 0.0]

x_initial = [0.8]

prob = ODEProblem(airfoils, x_initial, tspan, parameters)

sol = DifferentialEquations.solve(prob, reltol=1e-8)

answer = parsesolution(af, sol, parameters)

plot(answer[1,:], answer[2,:])