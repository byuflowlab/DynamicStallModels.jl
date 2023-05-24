using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations


dsm = DynamicStallModels


path = dirname(@__FILE__)
cd(path)


file = "../../polars/Larsen_Vertol_Polar.csv"
file_2 = "../../polars/Larsen_Vertol_Experimental.csv"
polar = readdlm(file, ',')
larsen_experimental = readdlm(file_2, ',')


c = 1.5
Vrel = 60

function uvec(t)
    return Vrel
end
function alphavec(t)
    c = 1.5
    Vrel = 60
    shift = 14.92
    amp = 4.85
    k = 0.062

    omega = 2*Vrel*k/c

    alf = shift + amp*sin(omega*t)
    
    return alf*(pi/180)
end

#############Functional Parameters###############
tspan = (0, 2.0)
xinit = [0.8, 0.8, 0.8, 0.8]
params = [uvec, 0.0, alphavec, 0.0, uvec, 0.0, alphavec, 0.0, uvec, 0.0, alphavec, 0.0, uvec, 0.0, alphavec, 0.0]
#################################################


#############Indicial Parameters################
tvec = range(0.0, 2.0, 1000)

aoa = alphavec.(tvec)
inflow_vec = ones(length(tvec)).*Vrel
################################################

functional_matrix = zeros(138, 12)
functional_time_matrix = zeros(138, 3)
indicial_matrix = zeros(1000, 12)



for i in 1:3
    ##################Functional Solve#################
    dsmodel_1 = Oye(Functional(), 1, i, 4.0)


    af_1_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.OSP())
    af_2_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.ADGSP())
    af_3_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.RSP())
    af_4_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.ADSP(1,1), eta=1.0)



    airfoils_1 = Array{Airfoil, 1}(undef , 4)
    airfoils_1[1] = af_1_Functional
    airfoils_1[2] = af_2_Functional
    airfoils_1[3] = af_3_Functional
    airfoils_1[4] = af_4_Functional



    prob = ODEProblem(airfoils_1, xinit, tspan, params)
    sol = DifferentialEquations.solve(prob, reltol=1e-8)

    answer = parsesolution(dsmodel_1, airfoils_1, sol, params)

    functional_matrix[:, 4(i-1)+1] = vcat(answer[2,:])
    functional_matrix[:, 4(i-1)+2] = vcat(answer[4,:])
    functional_matrix[:, 4(i-1)+3] = vcat(answer[6,:])
    functional_matrix[:, 4(i-1)+4] = vcat(answer[8,:])

    functional_time_matrix[:,i] = vcat(answer[1,:])
    #################################################

    ###################Indicial Solve#################

    dsmodel_2 = Oye(Indicial(), 1, i, 4.0)

    af_1_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.OSP())
    af_2_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.ADGSP())
    af_3_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.RSP())
    af_4_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.ADSP(1,1) , eta=1.0)

    airfoils_2 = Array{Airfoil, 1}(undef, 4)
    airfoils_2[1] = af_1_Indicial
    airfoils_2[2] = af_2_Indicial
    airfoils_2[3] = af_3_Indicial
    airfoils_2[4] = af_4_Indicial

    states, loads = solve_indicial(airfoils_2, tvec, inflow_vec, aoa)

    indicial_matrix[:, 4(i-1)+1] = loads[:, 1]
    indicial_matrix[:, 4(i-1)+2] = loads[:, 4]
    indicial_matrix[:, 4(i-1)+3] = loads[:, 7]
    indicial_matrix[:, 4(i-1)+4] = loads[:, 10]

end

#########Functional Plots###############
Hansen_Full_Sep = plot(functional_time_matrix[:,1].*180/pi, functional_matrix[:,1], xlabel = L"\alpha", ylabel = L"C_L", label = "OSP (Hansen)", linewidth =1.5, legend=:topleft, linestyle=:dash)
plot!(functional_time_matrix[:,1].*180/pi, functional_matrix[:,2], label = "ADGSP (Hansen)", linewidth = 1.5, linestyle=:dash)
plot!(functional_time_matrix[:,1].*180/pi, functional_matrix[:,3], label = "RSP (Hansen)", linewidth = 1.5, linestyle=:dash)
plot!(functional_time_matrix[:,1].*180/pi, functional_matrix[:,4], label = "ADSP (Hansen)", linewidth = 1.5, linestyle=:dash)
scatter!(larsen_experimental[:,1], larsen_experimental[:,2], label = "Experimental")


Hermite_Full_Sep = plot(functional_time_matrix[:,1].*180/pi, functional_matrix[:,5], xlabel = L"\alpha", ylabel = L"C_L", label = "OSP (Hermite)", linewidth =1.5, legend=:topleft)
plot!(functional_time_matrix[:,1].*180/pi, functional_matrix[:,6], label = "ADGSP (Hermite)", linewidth = 1.5, linestyle=:dash)
plot!(functional_time_matrix[:,1].*180/pi, functional_matrix[:,7], label = "RSP (Hermite)", linewidth = 1.5, linestyle=:dash)
plot!(functional_time_matrix[:,1].*180/pi, functional_matrix[:,8], label = "ADSP (Hermite)", linewidth = 1.5, linestyle=:dash)
scatter!(larsen_experimental[:,1], larsen_experimental[:,2], label = "Experimental")

Oye_Full_Sep = plot(functional_time_matrix[:,1].*180/pi, functional_matrix[:,9], xlabel = L"\alpha", ylabel = L"C_L", label = "OSP (Øye)", linewidth =1.5, legend=:topleft, linestyle=:dash)
plot!(functional_time_matrix[:,1].*180/pi, functional_matrix[:,10], label = "ADGSP (Øye)", linewidth = 1.5, linestyle=:dash)
plot!(functional_time_matrix[:,1].*180/pi, functional_matrix[:,11], label = "RSP (Øye)", linewidth = 1.5, linestyle=:dash)
plot!(functional_time_matrix[:,1].*180/pi, functional_matrix[:,12], label = "ADSP (Øye)", linewidth = 1.5, linestyle=:dash)
scatter!(larsen_experimental[:,1], larsen_experimental[:,2], label = "Experimental")
########################################

##################Indicial Plots################

Hansen_Full_Sep_2 = plot(aoa.*180/pi, indicial_matrix[:,1], xlabel = L"\alpha", ylabel = L"C_L", label = "OSP (Hansen)", linewidth =1.5, legend=:topleft)
plot!(aoa.*180/pi, indicial_matrix[:,2], label = "ADGSP (Hansen)", linewidth = 1.5)
plot!(aoa.*180/pi, indicial_matrix[:,3], label = "RSP (Hansen)", linewidth = 1.5, linestyle=:dash)
plot!(aoa.*180/pi, indicial_matrix[:,4], label = "ADSP (Hansen)", linewidth = 1.5)
scatter!(larsen_experimental[:,1], larsen_experimental[:,2], label = "Experimental")


Hermite_Full_Sep_2 = plot(aoa.*180/pi, indicial_matrix[:,5], xlabel = L"\alpha", ylabel = L"C_L", label = "OSP (Hermite)", linewidth =1.5, legend=:topleft)
plot!(aoa.*180/pi, indicial_matrix[:,6], label = "ADGSP (Hermite)", linewidth = 1.5, linestyle=:dash)
plot!(aoa.*180/pi, indicial_matrix[:,7], label = "RSP (Hermite)", linewidth = 1.5, linestyle=:dash)
plot!(aoa.*180/pi, indicial_matrix[:,8], label = "ADSP (Hermite)", linewidth = 1.5, linestyle=:dash)
scatter!(larsen_experimental[:,1], larsen_experimental[:,2], label = "Experimental")


Oye_Full_Sep_2 = plot(aoa.*180/pi, indicial_matrix[:,9], xlabel = L"\alpha", ylabel = L"C_L", label = "OSP (Øye)", linewidth =1.5, legend=:topleft, linestyle=:dash)
plot!(aoa.*180/pi, indicial_matrix[:,10], label = "ADGSP (Øye)", linewidth = 1.5, linestyle=:dash)
plot!(aoa.*180/pi, indicial_matrix[:,11], label = "RSP (Øye)", linewidth = 1.5, linestyle=:dash)
plot!(aoa.*180/pi, indicial_matrix[:,12], label = "ADSP (Øye)", linewidth = 1.5, linestyle=:dash)
scatter!(larsen_experimental[:,1], larsen_experimental[:,2], label = "Experimental")
################################################


