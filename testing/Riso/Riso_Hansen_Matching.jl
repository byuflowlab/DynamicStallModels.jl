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


function alphavec2(t)
    c = 0.55
    M = 0.1
    a = 343.0
    shift = 21.0
    amp = 4.0
    k = 0.1

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function alphavecdot2(t)
    c = 0.55
    M = 0.1
    a = 343.0
    shift = 21.0
    amp = 4.0
    k = 0.1

    v = M*a
    omega = k*2*v/c

    alf = amp*omega*cos(omega*t)
    return alf*(pi/180)
end

parameters2 = [Uvector, 0.0, alphavec2, alphavecdot2]
prob2 = ODEProblem(airfoils, x_initial, tspan, parameters2)
sol2 = DifferentialEquations.solve(prob2, reltol=1e-8)

answer2 = parsesolution(dsmodel, airfoils, sol2, parameters2)


function alphavec3(t)
    c = 0.55
    M = 0.1
    a = 343.0
    shift = 3.0
    amp = 4.0
    k = 0.1

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function alphavecdot3(t)
    c = 0.55
    M = 0.1
    a = 343.0
    shift = 3.0
    amp = 4.0
    k = 0.1

    v = M*a
    omega = k*2*v/c

    alf = amp*omega*cos(omega*t)
    return alf*(pi/180)
end

parameters3 = [Uvector, 0.0, alphavec3, alphavecdot3]
prob3 = ODEProblem(airfoils, x_initial, tspan, parameters3)
sol3 = DifferentialEquations.solve(prob3, reltol=1e-8)

answer3 = parsesolution(dsmodel, airfoils, sol3, parameters3)


plot(answer[1,50:end].*180/pi, answer[2,50:end], linewidth = 3, xlabel = "Angle of Attack (Degrees)", ylabel = "Cl", label = "DSM", legend=:bottomright,
 color=:blue, tickfontsize=10)

file2 = "../../polars/Hansen_6315_Blue_Full.csv"

Hansen = readdlm(file2, ',')

plot!(Hansen[:,1], Hansen[:,2], linewidth = 3, linestyle=:dash, label = "Hansen", color=:orange)

plot!(polar[55:90,1].*180/pi, polar[55:90,2], label = "Static", linestyle=:dashdot, linewidth =3)

file3 = "../../polars/Hansen_6315_Green.csv"

Hansen_Green = readdlm(file3, ',')

plot!(answer2[1,70:end].*180/pi, answer2[2,70:end], linewidth=3, label=:false, color=:blue)
plot!(Hansen_Green[:,1], Hansen_Green[:,2], linestyle=:dash, linewidth=3, label=:false, color=:orange)


file4 = "../../polars/Hansen_6315_Black.csv"
Hansen_Black = readdlm(file4, ',')

plot!(answer3[1,50:end].*180/pi, answer3[2,50:end], linewidth=3, label=:false, color=:blue)
plot!(Hansen_Black[:,1], Hansen_Black[:,2], linestyle=:dash, linewidth=3, label=:false, color=:orange)


#savefig("Riso_Hansen_Optimization_Grid")