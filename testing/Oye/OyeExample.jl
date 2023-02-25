#= OyeExample.jl
# Jacob Child
This is is a walk through example of how to use the Oye model, using the Vertol 23010-1.58 airfoil as an example.
All values are taken from the Larsen paper
=#

# Specify what packages we will need to use
using Revise, DynamicStallModels, DelimitedFiles, Plots, FLOWMath, OpenFASTsr, LaTeXStrings

#rename packages for convenience
dsm = DynamicStallModels
of = OpenFASTsr

# Get the path of the current file and change the working directory to that path
path = dirname(@__FILE__)
cd(path)

# Specify the airfoil and environmental variables, taken from pg 968
c = 1.5 # m, chord length
V = 60.0 #m/s, velocity 

# Extract Polar data from the file
VertolPolar = readdlm("C:/Users/child/Documents/Projects/FlowLab_DynamicStall/DynamicStallModels.jl/polars/extendedVertol_23010-1.58.dat", ',')
# Read in the airfoil data 
#Vertol = of.read_airfoilinput("../../data/airfoils/Vertol.dat") #read in the airfoil data using OpenFASTsr
#af = of.make_dsairfoil(Vertol) #make the airfoil into a DynamicStallModels airfoil
#af = dsm.simpleairfoil(VertolPolar) #? testing this versus the dsm.airfoil
af = dsm.airfoil(VertolPolar; A = 7.140)
airfoils = Array{Airfoil, 1}(undef, 1) #make an array of the type Airfoil struct
airfoils[1] = af #put the airfoil into the array

# Make the Oye model struct
dsmodel = Oye(Indicial(), 1, airfoils) #makes the struct, says it will solve it indicially and tha there is 1 airfoil

# Create time, velocity, and angle of attack vectors
tvec = range(0, 2.0, 1000) #time vector, these will be specific to the experimental data I am verifying against
Uvec = V.*ones(length(tvec)) #velocity vector across time

# Define the angle of attack function
function alpha(t)
    c = 1.5 #m
    v = 60.0 #m/s
    shift = 14.92 #degrees, this is the estimated mean angle of attack from pg 971, fig 9, plt d
    amp = 4.85 #degrees, amplitude of oscillation
    k = 0.062 #reduced frequency
    omega = k*2*v/c #rad/s, frequency of oscillation

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

# Create the angle of attack vector
alphavec = alpha.(tvec) #angle of attack vector across time, ie angle of attack at every time step

# Solve for the states and loads
states, loads = solve_indicial(dsmodel, [c], tvec, Uvec, alphavec) #solves the indicial model -> found under src/solve.jl

# Plot the states
stateplt = plot(tvec, states[:,1], leg=false, xaxis="Time (s)", yaxis="f")

# Plot the loads
cn = loads[:,1]
cn_static = af.cn.(alphavec)

cnplt = plot(xaxis="Time (s)", yaxis=L"C_n", leg=:topright)
plot!(tvec, cn, lab="DSM")
plot!(tvec, cn_static, lab="Static")

# Plot the whole cycle 
cyclecnplt = plot(xaxis="Angle of Attack (deg)", yaxis=L"C_n", leg=:topright)
plot!(alphavec*180/pi, cn, lab="DSM")
plot!(alphavec*180/pi, cn_static, lab="Static")

display(cyclecnplt)

