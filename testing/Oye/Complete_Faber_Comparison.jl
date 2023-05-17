using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations


dsm = DynamicStallModels


path = dirname(@__FILE__)
cd(path)


file = "../../polars/NACA_0015_Faber.csv"
polar = readdlm(file, ',')


c = 0.55
M = 0.11
a = 343.0
Vrel = M*a 


dsmodel = Oye(Functional(), 1, 3, 4.0)


af_1 = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.LSP())
af_2 = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.ADGSP())
af_3 = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.RSP())

airfoils = Array{Airfoil, 1}(undef, 3)
airfoils[1] = af_1
airfoils[2] = af_2
airfoils[3] = af_3


tspan = (0, 2.0)


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


parameters = [Uvector, 0.0, alpha, 0.0, Uvector, 0.0, alpha, 0.0, Uvector, 0.0, alpha, 0.0]
x_initial = [0.9, 0.9, 0.9]


prob = ODEProblem(airfoils, x_initial, tspan, parameters)
sol = DifferentialEquations.solve(prob, reltol=1e-8)


answer = parsesolution(dsmodel, airfoils, sol, parameters)


plot(answer[1,:], answer[2,:], xlabel=L"\mathrm{Angle ~of ~Attack ~(Radians)}", ylabel = L"C_L", label = "Larsen Separation Point", linewidth=1.5)
plot!(answer[3,:], answer[4,:], label = "AeroDyn Separation Point", linewidth=1.5)
plot!(answer[5,:], answer[6,:], label = "Riso Separation Point", linewidth=1.5)