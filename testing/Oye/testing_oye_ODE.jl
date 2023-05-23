using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations

dsm = DynamicStallModels


path = dirname(@__FILE__)
cd(path)


file = "../../polars/NACA_0015_Faber.csv" #this airfoil is from Faber's paper in the moderate stall section (NACA 0015)
polar_0015 = readdlm(file, ',')
polar = polar_0015


c = 0.55 #this block contains the chord length, mach number, and speed of sound
M = 0.11
a = 343.0
Vrel = M*a 


dsmodel = Oye(Functional(), 1, 2, 4.0) #this uses the Oye model with the functional solve. The flags inside, can be found in the Oye.jl file.
af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.OSP()) #the airfoil is created and will use the larsen/faber separation point function.


airfoils = Array{Airfoil, 1}(undef, 1) #a vector for the airfoils is created, and the single airfoil that we are evaluating is pushed in.
airfoils[1] = af


tspan = (0, 2.0) #this is the span of time values that the differential equations solver will evaluate over.


function Uvector(t) #these functions define the inflow velocity and angle of attack with respect to time for the airfoil
    return 0.11*343.0
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


parameters = [Uvector, 0.0, alpha, 0.0] #this is the parameter vector for the airfoil that will be used in the ODE solve
x_initial = [0.8] #this is the intial state value for the state equation (anything between 0 and 1 should be fine)


prob = ODEProblem(airfoils, x_initial, tspan, parameters)
sol = DifferentialEquations.solve(prob, reltol=1e-8)


answer= parsesolution(dsmodel, airfoils, sol, parameters) #parsesolution allows us to change the state values to the dynamic lift coefficients that we desire.
                                                          #parsesolution also gives the corresponding angle of attack values for the lift coefficients.


plot(answer[1,:].*180/pi, answer[2,:], xlabel = L"\mathrm{Angle~of~Attack~(Degrees)}", ylabel = L"C_L", label = "Oye", linewidth = 2)