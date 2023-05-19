using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations


dsm = DynamicStallModels 


path = dirname(@__FILE__)
cd(path)


file = "../../polars/NACA_0015_Faber.csv" #this polar is from Faber's paper from his moderate stall section. It is the NACA 0015 airfoil.
polar = readdlm(file, ',')


c = 0.55 #this is the chord length of the airfoil
M = 0.11 #this is the mach number of the airofil
a = 343.0 #this is the speed of sound
Vrel = M*a 


############# Functional Approach #############
dsmodel_1 = Oye(Functional(), 1, 2, 4.0) #this generates the the model that will be used for the functional solve. It uses coefficient of lift,
                                         #the larsen/faber method of finding fully separated lift, and states that the time constant is 4.0


af_1_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.LSP()) #this section creates the airfoil structs for three different approaches to separation point.
af_2_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.ADGSP())
af_3_Functional = dsm.make_airfoil(polar, dsmodel_1, c; sfun=dsm.RSP())


airfoils_1 = Array{Airfoil, 1}(undef, 3) #this block of code creates the airfoil vector and pushes the airfoil structs into this vector
airfoils_1[1] = af_1_Functional
airfoils_1[2] = af_2_Functional
airfoils_1[3] = af_3_Functional


tspan = (0, 2.0) #this creates the span of time that the differential equation solver will evaluate over


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


parameters = [Uvector, 0.0, alpha, 0.0, Uvector, 0.0, alpha, 0.0, Uvector, 0.0, alpha, 0.0] #this is the parameter vector that will be passed into the solver
x_initial = [0.9, 0.9, 0.9] #this is the initial state value for the state rate equation (anything between 0 and 1 should be fine)


prob = ODEProblem(airfoils_1, x_initial, tspan, parameters) 
sol = DifferentialEquations.solve(prob, reltol=1e-8)


answer = parsesolution(dsmodel_1, airfoils_1, sol, parameters) #the parsesolution function turns the state values into the dynamic lift coefficients and
                                                               #also gives the corresponding angle of attacks 


#plot(answer[1,:].*180/pi, answer[2,:], xlabel=L"\mathrm{Angle ~of ~Attack ~(Degrees)}", ylabel = L"C_L", label = "Larsen Separation Point", linewidth=1.5)
#plot!(answer[3,:].*180/pi, answer[4,:], label = "AeroDyn Separation Point", linewidth=1.5)
#Functional_Plot = plot!(answer[5,:].*180/pi, answer[6,:], label = "Riso Separation Point", linewidth=1.5)


########### Indicial Approach ############
dsmodel_2 = Oye(Indicial(), 1, 2, 4.0) #this creates the model that will be used for this indicial solve. It has the same conditions as the functional solve.


af_1_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.LSP()) #this block is the same set up as the functional solve.
af_2_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.ADGSP())
af_3_Indicial = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.RSP())


airfoils_2 = Array{Airfoil, 1}(undef, 3) #the airfoils for the indicial solve are pushed into this airofil vector in this block
airfoils_2[1] = af_1_Indicial
airfoils_2[2] = af_2_Indicial
airfoils_2[3] = af_3_Indicial


tvec = range(0.0, 1.0, 1000) #this range gives the time interval that will be evaluated over and gives the size of the time steps


aoa = alpha.(tvec) #using the angle of attack function from earlier, a vector of the angle of attacks is created using the time range
Uvec = Vrel.*ones(length(tvec)) #like the aoa generation, this inflow vector is created using the time range.


states , loads = solve_indicial(airfoils_2, tvec, Uvec, aoa) 


plot(aoa, loads[:,1], xlabel = L"\mathrm{Angle~of~Attack~(Radians)}", ylabel = L"C_L", label = "Larsen Separation Point", linewidth =1.5)
plot!(aoa, loads[:,4], label =  "AeroDyn Separation Point",  linewidth =1.5)
plot!(aoa, loads[:,7], label = "Riso Separation Point",  linewidth =1.5)
