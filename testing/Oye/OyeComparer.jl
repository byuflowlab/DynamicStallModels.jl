#OyeComparer.jl
#To compare DSM Oye the modeled data to both the Larsen Oye model
# and the Larsen experimental data

# The code below is copied from OyeExample.jl to get the data
    # later it can be all be replaced with a file location for the data,
    # but for now things might still need to be changed

    # Specify what packages we will need to use
println("OyeComparer is running!")
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
af = dsm.airfoil(VertolPolar; A = .07, sfun=dsm.LSP()) #A= 7.14
af = dsm.update_airfoil(af; alphasep=[af.alphasep[1], 32.0*pi/180] ) #update the airfoil with Larsen's separation point
airfoils = Array{Airfoil, 1}(undef, 1) #make an array of the type Airfoil struct
airfoils[1] = af #put the airfoil into the array

# Make the Oye model struct
#cflag::Int - A flag to apply the separation delay to the coefficient of 1) lift, 2) normal force. 
#version::Int - A flag to say whether to use 1) Hansen 2008, or 2) Faber 2018's implementation of the model, or 3) BeddoesLeishman, or 4) Larsen's 2007
dsmodel = Oye(Indicial(), 1, airfoils,1,4) #makes the struct, says it will solve it indicially and that there is 1 airfoil
#! if I call the above with a 3 (technically beddoesLeishman) it errors in Oye 51 or solve 56, or solve 11?
# Create time, velocity, and angle of attack vectors
tvec = range(0, 2.0, 1000) #time vector, these will be specific to the experimental data I am verifying against
Uvec = V.*ones(length(tvec)) #velocity vector across time

# Define the angle of attack function
function alpha(t)
    c = 1.5 #m
    v = 60.0 #m/s
    shift = 15.0*pi/180 #14.92 #degrees, this is the estimated mean angle of attack from pg 971, fig 9, plt d
    amp = 4.85*pi/180 #degrees, amplitude of oscillation
    k = .062 #0.062 #reduced frequency
    omega = k*2*v/c #rad/s, frequency of oscillation

    alf = shift + amp*sin(omega*t)
    return alf
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
#I am going to use a bunch of vectors rather than matrices as Julia has good vector functions
# Import the data #! this is hardcoded with my info for now!
LarsenOyeDyn = readdlm("C:\\Users\\child\\Documents\\Projects\\FlowLab_DynamicStall\\DynamicStallModels.jl\\data\\Larsen2007\\Oye\\WPDLarsenFig9dDataTrimmed.csv", ',') #dynamic oye data
LarsenStatic = readdlm("C:\\Users\\child\\Documents\\Projects\\FlowLab_DynamicStall\\DynamicStallModels.jl\\data\\Larsen2007\\Oye\\WPDLarsenFig9DataStatic.csv", ',') #static data
LarsenOyeDynDeg = vec(LarsenOyeDyn[:,1]) #angle of attack in degrees
LarsenOyeDynCn = LarsenOyeDyn[:,2] #normal force coefficient
LarsenStaticDeg = LarsenStatic[:,1] #angle of attack in degrees
LarsenStaticCn = LarsenStatic[:,2] #normal force coefficient
#TODO: extract (web plot digitizer) larsen experimental data LarsenOyeExp = readdlm("C:\\Users\\child\\Documents\\Projects\\FlowLab_DynamicStall\\DynamicStallModels.jl\\data\\Larsen2007\\Oye\\WPDLarsenFig9DataExp.csv", ',')

#Split the dynamic data into the upstroke (while the angle is increasing) and downstroke (while the angle is decreasing)
#println("length of LarsenOyeDyn: ", length(LarsenOyeDyn))
LarsenOyeDynUpDeg = vec(zeros(length(LarsenOyeDyn[:,1]), 1))
LarsenOyeDynUpCn = vec(zeros(length(LarsenOyeDyn[:,1]), 1))
LarsenOyeDynDownDeg = vec(zeros(length(LarsenOyeDyn[:,1]), 1))
LarsenOyeDynDownCn = vec(zeros(length(LarsenOyeDyn[:,1]), 1))
# for loop to split and extract the data 
LarsenOyeDynUpDeg[1] = LarsenOyeDynDeg[1] #!something doesn't work in this loop!
LarsenOyeDynUpCn[1] = LarsenOyeDynCn[1]
for i in 2:length(LarsenOyeDynDeg)
    #println(i)
    if LarsenOyeDynDeg[i] > LarsenOyeDynDeg[i-1] #? what is the blue warning? each index or axes?
        LarsenOyeDynUpDeg[i] = LarsenOyeDyn[i]
        LarsenOyeDynUpCn[i] = LarsenOyeDynCn[i]
        #println("if statement ", LarsenOyeDynUp)
    elseif LarsenOyeDyn[i] < LarsenOyeDyn[i-1]
        LarsenOyeDynDownDeg[i] = LarsenOyeDynDeg[i]
        LarsenOyeDynDownCn[i] = LarsenOyeDynCn[i]
        #println("elseif ", LarsenOyeDynDown)
    else 
        #println("I am here")
    end
end
#Clean the zeros from the vectors
deleteat!(LarsenOyeDynUpDeg, findall(x->x==0, LarsenOyeDynUpDeg))
deleteat!(LarsenOyeDynUpCn, findall(x->x==0, LarsenOyeDynUpCn))
deleteat!(LarsenOyeDynDownDeg, findall(x->x==0, LarsenOyeDynDownDeg))
deleteat!(LarsenOyeDynDownCn, findall(x->x==0, LarsenOyeDynDownCn))

#plot to check
scatter(alphavec*180/pi, cn, label = "DSM")
scatter!(LarsenOyeDynUpDeg, LarsenOyeDynUpCn, label = "Larsen Downstroke")
scatter!(LarsenOyeDynDownDeg, LarsenOyeDynDownCn, label = "Larsen Upstroke", legend=:topleft)