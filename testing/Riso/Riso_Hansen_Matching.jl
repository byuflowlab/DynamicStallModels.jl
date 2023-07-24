using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations

dsm = DynamicStallModels


path = dirname(@__FILE__)
cd(path)

file = "../../polars/Hansen_6315.csv"

polar = readdlm(file , ',')

c = 0.55
M = 0.1
v = 343
Vrel = M*v

dsmodel = Riso(Functional(), [0.27469922114926987, 0.05433173899834567] , [0.11224679387440323, 0.45235838971138226], [0.8173872939878122, 3.671051417177959]) 
af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.RSP())

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af

tspan = (0.0, 2.0)

function Uvector(t)
    return 0.1*343
end

function alphavec(t)
    c = 0.55
    M = 0.1
    a = 343.0
    shift = 12.0
    amp = 4.0
    k = 0.1

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function alphavecdot(t)
    c = 0.55
    M = 0.1
    a = 343.0
    shift = 12.0
    amp = 4.0
    k = 0.1

    v = M*a
    omega = k*2*v/c

    alf = amp*omega*cos(omega*t)
    return alf*(pi/180)
end

parameters = [Uvector, 0.0, alphavec, alphavecdot]
x_initial = [0.0, 0.0, 0.0, 0.0]

prob = ODEProblem(airfoils, x_initial, tspan, parameters)
sol = DifferentialEquations.solve(prob, reltol=1e-8)

answer = parsesolution(dsmodel, airfoils, sol, parameters)

plot(answer[1,:].*180/pi, answer[2,:], linewidth = 2, xlabel = "Angle of Attack (Degrees)", ylabel = "Cl", label = "Optimized DSM")

file2 = "../../polars/Hansen_6315_Blue_Full.csv"

Hansen = readdlm(file2, ',')

plot!(Hansen[:,1], Hansen[:,2], linewidth = 2, linestyle=:dash, label = "Hansen's Results")

plot!(polar[55:80,1].*180/pi, polar[55:80,2], label = "Static")

#savefig("Riso_Hansen_Optimization")