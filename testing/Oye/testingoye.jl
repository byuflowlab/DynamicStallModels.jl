using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, OpenFASTsr, LaTeXStrings

dsm = DynamicStallModels #renames the package to make using it easier
of = OpenFASTsr

path = dirname(@__FILE__) #gets the path of the current file
cd(path) #changes the working directory to the path of the current file

c = 0.1 #chord 

M = 0.379   #Mach number
a = 343.0   #speed of sound
Vrel = M*a #60

# polar = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/learning/exploring/exampledata/NACA4412.dat", '\t'; skipstart=3) 
# af = airfoil(polar) #Todo: This constructor is broken.

du21_a17 = of.read_airfoilinput("../../data/airfoils/DU40_A17.dat") 
af = of.make_dsairfoil(du21_a17) #Todo: I think this polar might be my problem. I should try a different polar.... which means that I need to fix the constructor. :| 

# af = update_airfoil(af, A=[4.0], dcndalpha=6.320368333107256, alpha0=-0.0033903071711640564)

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af


dsmodel = Oye(Indicial(), 1, airfoils) #This makes the dsmodel Oye struct, it is a method on a struct and hasn't run the function yet


tvec = range(0, 0.05, 100) #0:0.001:0.05 #time vector, these will be specific to the experimental data I am verifying against
Uvec = Vrel.*ones(length(tvec)) #velocity vector across time

function alpha(t)
    c = 0.1
    M = 0.379
    a = 343.0
    shift = 10.3
    amp = 8.1
    k = 0.075 #reduced frequency 

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

alphavec = alpha.(tvec) #angle of attack vector across time, ie angle of attack at every time step

states, loads = solve_indicial(dsmodel, [c], tvec, Uvec, alphavec) #solves the indicial model-> found under src/solve.jl


stateplt = plot(tvec, states[:,1], leg=false, xaxis="time (s)", yaxis="f")
# display(stateplt)


cn = loads[:,1]
cn_static = af.cn.(alphavec)

cnplt = plot( xaxis="time (s)", yaxis=L"C_n", leg=:topright)
plot!(tvec, cn, lab="DSM")
plot!(tvec, cn_static, lab="Static")
# display(cnplt) 


cyclecnplt = plot(xaxis="Angle of Attack (deg)", yaxis=L"C_n", leg=:topright)
plot!(alphavec.*(180/pi), cn, lab="DSM")
plot!(alphavec.*(180/pi), cn_static, lab="Static")
display(cyclecnplt) 

#=
I think that I have the models matching, the difference is that I'm using Cn in this model and Cl in the other. 
=#


nothing