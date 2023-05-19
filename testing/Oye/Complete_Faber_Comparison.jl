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

############# Functional Approach #############
dsmodel_1 = Oye(Functional(), 1, 2, 4.0)


af_1_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.LSP())
af_2_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.ADGSP())
af_3_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.RSP())

airfoils_1 = Array{Airfoil, 1}(undef, 3)
airfoils_1[1] = af_1_Functional
airfoils_1[2] = af_2_Functional
airfoils_1[3] = af_3_Functional


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


prob = ODEProblem(airfoils_1, x_initial, tspan, parameters)
sol = DifferentialEquations.solve(prob, reltol=1e-8)


answer = parsesolution(dsmodel_1, airfoils_1, sol, parameters)


#plot(answer[1,:].*180/pi, answer[2,:], xlabel=L"\mathrm{Angle ~of ~Attack ~(Degrees)}", ylabel = L"C_L", label = "Larsen Separation Point", linewidth=1.5)
#plot!(answer[3,:].*180/pi, answer[4,:], label = "AeroDyn Separation Point", linewidth=1.5)
#Functional_Plot = plot!(answer[5,:].*180/pi, answer[6,:], label = "Riso Separation Point", linewidth=1.5)


########### Indicial Approach ############

dsmodel_2 = Oye(Indicial(), 1, 2, 4.0)


af_1_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.LSP())
af_2_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.ADGSP())
af_3_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.RSP())

airfoils_2 = Array{Airfoil, 1}(undef, 3)
airfoils_2[1] = af_1_Indicial
airfoils_2[2] = af_2_Indicial
airfoils_2[3] = af_3_Indicial

tvec = range(0.0, 1.0, 1000)

aoa = alpha.(tvec)

Uvec = Vrel.*ones(length(tvec))

states , loads = solve_indicial(airfoils_2, tvec, Uvec, aoa)

plot(aoa, loads[:,1], xlabel = L"\mathrm{Angle~of~Attack~(Radians)}", ylabel = L"C_L", label = "Larsen Separation Point", linewidth =1.5)
plot!(aoa, loads[:,4], label =  "AeroDyn Separation Point",  linewidth =1.5)
plot!(aoa, loads[:,7], label = "Riso Separation Point",  linewidth =1.5)
