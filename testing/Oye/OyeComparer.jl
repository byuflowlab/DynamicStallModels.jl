#OyeComparer.jl
#To compare DSM Oye the modeled data to both the Larsen Oye model
# and the Larsen experimental data

# The code below is copied from OyeExample.jl to get the data
    # later it can be all be replaced with a file location for the data,
    # but for now things might still need to be changed

    # Specify what packages we will need to use
using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, OpenFASTsr, LaTeXStrings

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

#TODO plots are commented out for now
cnplt = plot(xaxis="Time (s)", yaxis=L"C_n", leg=:topright)
#plot!(tvec, cn, lab="DSM")
#plot!(tvec, cn_static, lab="Static")

# Plot the whole cycle 
cyclecnplt = plot(xaxis="Angle of Attack (deg)", yaxis=L"C_n", leg=:topright)
#plot!(alphavec*180/pi, cn, lab="DSM")
#plot!(alphavec*180/pi, cn_static, lab="Static")

#display(cyclecnplt)

# Now we will compare the results to the Larsen Oye model
# Import the data #! this is hardcoded with my info for now!
LarsenOyeDyn = readdlm("C:\\Users\\child\\Documents\\Projects\\FlowLab_DynamicStall\\DynamicStallModels.jl\\data\\Larsen2007\\Oye\\WPDLarsenFig9dDataTrimmed.csv", ',') #dynamic oye data
LarsenStatic = readdlm("C:\\Users\\child\\Documents\\Projects\\FlowLab_DynamicStall\\DynamicStallModels.jl\\data\\Larsen2007\\Oye\\WPDLarsenFig9DataStatic.csv", ',') #static data
#TODO: extract (web plot digitizer) larsen experimental data LarsenOyeExp = readdlm("C:\\Users\\child\\Documents\\Projects\\FlowLab_DynamicStall\\DynamicStallModels.jl\\data\\Larsen2007\\Oye\\WPDLarsenFig9DataExp.csv", ',')

#Split the dynamic data into the upstroke (while the angle is increasing) and downstroke (while the angle is decreasing)
#println("length of LarsenOyeDyn: ", length(LarsenOyeDyn))
LarsenOyeDynUp = zeros(length(LarsenOyeDyn[:,1]), 2)
LarsenOyeDynDown = zeros(length(LarsenOyeDyn[:,1]), 2)
# for loop to split and extract the data 
LarsenOyeDynUp[1,1] = LarsenOyeDyn[1,1]
for i in 2:length(LarsenOyeDyn[:,1])
    println(i)
    if LarsenOyeDyn[i,1] > LarsenOyeDyn[i-1,1] #? what is the blue warning? each index or axes?
        LarsenOyeDynUp[i,1] = LarsenOyeDyn[i,1]
        LarsenOyeDynUp[i,2] = LarsenOyeDyn[i,2]
    else
        LarsenOyeDynDown[i,1] = LarsenOyeDyn[i,1]
        LarsenOyeDynDown[i,2] = LarsenOyeDyn[i,2]
    end
end