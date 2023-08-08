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

dsmodel = Riso(Functional(), [0.11901046494134122, 0.2964635220443628] , [0.03833505254532422, 1.0129092704149298], [1.2267722555383287, 6.257533490189598]) 
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
    amp = 10.64
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
    amp = 10.64
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

answer = parsesolution(dsmodel, airfoils, sol, parameters)

plot(polar[52:85,1].*180/pi, polar[52:85,2], linestyle=:dashdot, linewidth=3, label = "Static (NACA 0030)", color=:black)

plot!(answer[1,:].*180/pi, answer[2,:], xlabel = "Angle of Attack (Degrees)", ylabel = "Cl", linewidth=3, label = "DSM", color=:blue)

file2 = "../../polars/Faber_Riso_NACA_0030.csv"
Faber_Riso = readdlm(file2, ',')

plot!(Faber_Riso[:,1], Faber_Riso[:,2], linestyle=:dash, linewidth=3, label="Faber", color=:orange)
savefig("Faber_Riso_Matching_0030")