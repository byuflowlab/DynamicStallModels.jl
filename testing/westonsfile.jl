using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations
dsm = DynamicStallModels

path = dirname(@__FILE__)
cd(path)

# file = "../../polars/Larsen_Vertol_23010.csv"
file = "/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/data/polars/extendedVertol 23010-1.58.dat"
polar = readdlm(file, ',')
polar[:,1] = polar[:,1].*pi/180
# matrix = zeros(32, 4)
# matrix[1:32, 1] = polar[:,1]
# matrix[1:32, 2] = polar[:,2]
n, m = size(polar)
polar = hcat(polar, zeros(n))
c = 1.5
V = 60

# dsmodel = Larsen(Functional(), [3.64, 24, 8, 6] , [0.165, 0.335], 14.75*pi/180)
dsmodel = Oye(Functional(), 1, 1, 4)
af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.RSP())
af_new = dsm.update_airfoil(af;  dcldalpha=6.71329364062, alpha0=0.0154469632272)
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af_new
tspan = (0.0, 2.0)

error("stop!")

function Uvector(t)
    return 60
end

function alphavec(t)
    c = 1.5
    M = 0.1749
    a = 343.0
    shift = 14.92
    amp = 4.85
    k = 0.062
    v = M*a
    omega = k*2*v/c
    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function alphavecdot(t)
    c = 1.5
    M = 0.1749
    a = 343.0
    shift = 14.92
    amp = 4.85
    k = 0.062
    v = M*a
    omega = k*2*v/c
    alf = amp*omega*cos(omega*t)
    return alf*(pi/180)
end

# parameters = [Uvector, 0.0, alphavec, alphavecdot]
# x_initial = [0.0, 0.0, 0.0, 0.0]
# prob = ODEProblem(airfoils, x_initial, tspan, parameters)
# sol = DifferentialEquations.solve(prob, reltol=1e-8)
# new = Array(sol)
# answer = parsesolution(dsmodel, airfoils, sol, parameters)
# plot(answer[1,:].*180/pi, answer[2,:])
# file2 = "../../polars/Larsen_Vertol_23010_Comparison.csv"
# compare = readdlm(file2, ',')
# plot!(compare[:,1], compare[:,2])
