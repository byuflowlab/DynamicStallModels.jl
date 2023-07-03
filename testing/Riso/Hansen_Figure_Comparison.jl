using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations

dsm = DynamicStallModels


path = dirname(@__FILE__)
cd(path)

file = "../../polars/Hansen_Figure3.csv"

polar = readdlm(file , ',')

c = 0.55
M = 0.11
v = 343
Vrel = M*v

dsmodel = Riso(Functional(), [0.294, 0.331] , [0.0664, 0.3266], [1.5, 6]) 
af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.RSP())

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af

tspan = (0.0, 0.224)

function Uvector(t)
    return 0.11*343
end

function alphavec(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = -3.0
    amp = 46.0
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
    shift = -3.0
    amp = 46.0
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

plot(answer[1,:].*180/pi, stuff[1,:], linewidth = 2, label = "Static", xlabel = "Angle of Attack (Degrees)", legend =:right)
plot!(answer[1,:].*180/pi, more[1,:], linewidth = 2, label = "Separation Point (DSM)")
plot!(answer[1,:].*180/pi, extra[1,:], linewidth=2, label = "Fully Separated Lift (DSM)")

file2 = "../../polars/Hansen_Fig3_FullSep.csv"
Hansen_Full_Sep = readdlm(file2, ',')

plot!(Hansen_Full_Sep[:,1], Hansen_Full_Sep[:,2], linewidth =2, linestyle=:dash, label = "Fully Separated Lift (Hansen)")

file3 = "../../polars/Hansen_Fig3_SepPoint.csv"
Hansen_SepPoint = readdlm(file3, ',')

plot!(Hansen_SepPoint[:,1], Hansen_SepPoint[:,2], linewidth =2, linestyle=:dash, label = "Separation Point (Hansen)")

savefig("Hansen_Figure_Comparison")