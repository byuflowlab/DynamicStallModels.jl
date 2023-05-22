using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations


dsm = DynamicStallModels 


path = dirname(@__FILE__)
cd(path)


file = "../../polars/NACA_0015_Faber.csv" #this polar is from Faber's paper from his moderate stall section. It is the NACA 0015 airfoil.
file_2 = "../../polars/Faber_0015_Oye_Top_Half.csv"
file_3 = "../../polars/Faber_0015_Oye_Bottom_Half.csv"
polar = readdlm(file, ',')
Faber_Top = readdlm(file_2, ',')
Faber_Bottom = readdlm(file_3, ',')


c = 0.55 #this is the chord length of the airfoil
M = 0.11 #this is the mach number of the airofil
a = 343.0 #this is the speed of sound
Vrel = M*a 



function Uvector(t) #these two functions define what the inflow velocity and angle of attack are with respect to time for the airfoils
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




######Indicial Parameters##########
tvec = range(0.0, 1.0, 1000) #this range gives the time interval that will be evaluated over and gives the size of the time steps


aoa = alpha.(tvec) #using the angle of attack function from earlier, a vector of the angle of attacks is created using the time range
Uvec = Vrel.*ones(length(tvec)) #like the aoa generation, this inflow vector is created using the time range.
###################################




########Functional Parameters########
tspan = (0, 2.0) #this creates the span of time that the differential equation solver will evaluate over


parameters = [Uvector, 0.0, alpha, 0.0, Uvector, 0.0, alpha, 0.0, Uvector, 0.0, alpha, 0.0, Uvector, 0.0, alpha, 0.0] #this is the parameter vector that will be passed into the solver
x_initial = [0.9, 0.9, 0.9, 0.9] #this is the initial state value for the state rate equation (anything between 0 and 1 should be fine)
#####################################




Indicial_Matrix = zeros(1000, 12) 
Functional_Matrix = zeros(377, 12)
Functional_Time_Matrix = zeros(377, 3)

for i in 1:3
    ############# Functional Approach #############
    dsmodel_1 = Oye(Functional(), 1, i, 4.0) #this generates the the model that will be used for the functional solve. It uses coefficient of lift,
                                            #the larsen/faber method of finding fully separated lift, and states that the time constant is 4.0


    af_1_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.LSP()) #this section creates the airfoil structs for three different approaches to separation point.
    af_2_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.ADGSP())
    af_3_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.RSP())
    af_4_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.ADSP(1,1), eta=1.0)
    #af_5_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.BLSP())


    airfoils_1 = Array{Airfoil, 1}(undef, 4) #this block of code creates the airfoil vector and pushes the airfoil structs into this vector
    airfoils_1[1] = af_1_Functional
    airfoils_1[2] = af_2_Functional
    airfoils_1[3] = af_3_Functional
    airfoils_1[4] = af_4_Functional
    #airfoils_1[5] = af_5_Functional


    prob = ODEProblem(airfoils_1, x_initial, tspan, parameters) 
    sol = DifferentialEquations.solve(prob, reltol=1e-8)

    


    answer = parsesolution(dsmodel_1, airfoils_1, sol, parameters) #the parsesolution function turns the state values into the dynamic lift coefficients and
                                                                   #also gives the corresponding angle of attacks 


    Functional_Matrix[:, 4(i-1)+1] = vcat(answer[2, :])
    Functional_Matrix[:, 4(i-1)+2] = vcat(answer[4, :])
    Functional_Matrix[:, 4(i-1)+3] = vcat(answer[6, :])
    Functional_Matrix[:, 4(i-1)+4] = vcat(answer[8, :])

    Functional_Time_Matrix[:,i] = vcat(answer[1,:])

    #plot(answer[1,:].*180/pi, answer[2,:], xlabel=L"\mathrm{Angle ~of ~Attack ~(Degrees)}", ylabel = L"C_L", label = "Larsen Separation Point", linewidth=1.5)
    #plot!(answer[3,:].*180/pi, answer[4,:], label = "AeroDyn Separation Point", linewidth=1.5)
    #Functional_Plot = plot!(answer[5,:].*180/pi, answer[6,:], label = "Riso Separation Point", linewidth=1.5)


    ########### Indicial Approach ############
    dsmodel_2 = Oye(Indicial(), 1, i, 4.0) #this creates the model that will be used for this indicial solve. It has the same conditions as the functional solve.


    af_1_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.LSP()) #this block is the same set up as the functional solve.
    af_2_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.ADGSP())
    af_3_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.RSP())
    af_4_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.ADSP(1,1), eta=1.0)
    #af_5_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.BLSP())

    airfoils_2 = Array{Airfoil, 1}(undef, 4) #the airfoils for the indicial solve are pushed into this airofil vector in this block
    airfoils_2[1] = af_1_Indicial
    airfoils_2[2] = af_2_Indicial
    airfoils_2[3] = af_3_Indicial
    airfoils_2[4] = af_4_Indicial
    #airfoils_2[5] = af_5_Indicial


    states , loads = solve_indicial(airfoils_2, tvec, Uvec, aoa) 

    Indicial_Matrix[:, 4(i-1)+1] = loads[:,1]
    Indicial_Matrix[:, 4(i-1)+2] = loads[:,4]
    Indicial_Matrix[:, 4(i-1)+3] = loads[:,7]
    Indicial_Matrix[:, 4(i-1)+4] = loads[:,10]


    #plot(aoa, loads[:,1], xlabel = L"\mathrm{Angle~of~Attack~(Radians)}", ylabel = L"C_L", label = "Larsen Separation Point", linewidth =1.5)
    #plot!(aoa, loads[:,4], label =  "AeroDyn Separation Point",  linewidth =1.5)
    #plot!(aoa, loads[:,7], label = "Riso Separation Point",  linewidth =1.5)
    #plot!(aoa, loads[:,10], label = "ADSP",  linewidth =1.5)
    #plot!(aoa, loads[:,13], label = "BLSP",  linewidth =1.5)
end

############Indicial Plots##############


#plot(aoa.*180/pi , Indicial_Matrix[: , 1] , xlabel = L"\alpha" , ylabel = L"C_L", label = "LSP (Hansen)", linewidth = 1.5, legend=:topleft)
#plot!(aoa , Indicial_Matrix[: , 2] , label = "ADGSP (Hansen Fully Separated)")
#plot!(aoa.*180/pi , Indicial_Matrix[: , 3] , label = "RSP (Hansen)", linewidth = 1.5)
#scatter!(Faber_Top[:,1], Faber_Top[:,2], color=:black, label = "Faber's Øye Results")
#scatter!(Faber_Bottom[:,1], Faber_Bottom[:,2], color=:black, label = false)
#Indicial_Plot_Oye_Hansen = plot!(aoa.*180/pi , Indicial_Matrix[: , 4] , label = "ADSP (Hansen)", linewidth = 1.5)
#plot!(aoa.*180/pi , Indicial_Matrix[: , 5] , label = "LSP (Hermite Fully Separated)", linewidth = 1.5)
#plot!(aoa.*180/pi , Indicial_Matrix[: , 6] , label = "ADGSP (Hermite)", linewidth = 1.5)
#plot!(aoa.*180/pi , Indicial_Matrix[: , 7] , label = "RSP (Hermite)", linewidth = 1.5)
#plot!(aoa.*180/pi, Indicial_Matrix[: , 8] , label = "ADSP (Hermite)", linewidth = 1.5)
#plot!(aoa.*180/pi , Indicial_Matrix[: , 9] , label = "LSP (Oye Fully Separated)", linewidth = 1.5)
#plot!(aoa.*180/pi , Indicial_Matrix[: , 10] , label = "ADGSP (Øye)", linewidth = 1.5)
#plot!(aoa.*180/pi , Indicial_Matrix[: , 11] , label = "RSP (Øye)", linewidth = 1.5)
#plot!(aoa.*180/pi , Indicial_Matrix[: , 12] , label = "ADSP (Øye)", linewidth = 1.5)

#savefig("Indicial_Oye_Hansen_Faber_ADGSP")

##########################################

plot(Functional_Time_Matrix[:,1].*180/pi, Functional_Matrix[:, 9],  xlabel = L"\alpha" , ylabel = L"C_L", label = "LSP (Øye)", linewidth = 1.5, legend=:topleft)
plot!(Functional_Time_Matrix[:,1].*180/pi, Functional_Matrix[:, 10], label = "ADGSP (Øye)", linewidth = 1.5)
plot!(Functional_Time_Matrix[:,1].*180/pi, Functional_Matrix[:, 11], label = "RSP (Øye)", linewidth = 1.5)
plot!(Functional_Time_Matrix[:,1].*180/pi, Functional_Matrix[:, 12], label = "ADSP (Øye)", linewidth = 1.5)
scatter!(Faber_Top[:,1], Faber_Top[:,2], color=:black, label = "Faber's Øye Results")
scatter!(Faber_Bottom[:,1], Faber_Bottom[:,2], color=:black, label = false)

savefig("Functional_Oye_Oye_Faber")
