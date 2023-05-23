using DynamicStallModels, DelimitedFiles


dsm = DynamicStallModels


file = "../../polars/NACA_0015_Faber.csv" #the airfoil used is from Faber's paper in the moderate stall section (NACA 0015)
polar_0015 = readdlm(file , ',')


dsmodel = Oye(Indicial(), 1, 2, 4.0) #this generates the model that will be used. We will be using the Oye method with the indical solve approach.
                                     #the description of the flags found inside can be found in the Oye.jl file.


c = 0.55 #this block gives the chord length, mach number, and the speed of sound.
M = 0.11
a = 343.0
Vrel = M*a 


af = dsm.make_airfoil(polar_0015, dsmodel, c; sfun=dsm.OSP()) #the airfoil struct is created here and specifically calls for the larsen/faber separation point method.
airfoils = Array{Airfoil, 1}(undef, 1) #the cevtor for the airfoils is created, and the single airfoil struct is pushed into this vector.
airfoils[1] = af


tvec = range(0, 1.0, 1000) #this gives the range of time values that will be evaluated over and also gives the magnitude of the time step.


function alpha(t) #this function gives the angle of attack of the airfoil as a function of time.
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


alphavec = alpha.(tvec) #using the angle of attack function, a vector of angle of attack values is created using the time range.
Uvec = Vrel.*ones(length(tvec)) #a similar vector to the angle of attack vector is created but for the inflow velocity.


states, loads = solve_indicial(airfoils, tvec, Uvec, alphavec)


cn = loads[:,1] #this block can either change the calculated coefficient values to be with respect to lift or the normal force.
if dsmodel.cflag == 2 #with the changes that we have made, though, I'm not sure if this block is necessary anymore.
    cn_static = af.cn.(alphavec)
else
    cn_static = af.cl.(alphavec)
end


using Plots, LaTeXStrings
stateplt = plot(tvec, states[:,1], leg=false, xaxis="time (s)", yaxis="f")
cyclecnplt = plot(xaxis="Angle of Attack (deg)", yaxis=L"C_L", leg=:bottomright)
plot!(alphavec.*(180/pi), cn, lab="DSM")
plot!(alphavec.*(180/pi), cn_static, lab="Static")
display(cyclecnplt) 