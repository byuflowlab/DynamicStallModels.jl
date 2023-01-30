#= OyeExample.jl
# Jacob Child
This is is a walk through example of how to use the Oye model, using the Vertol 23010-1.58 airfoil as an example.
All values are taken from the Larsen paper
=#

# Specify what packages we will need to use
using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, OpenFastsr, LaTeXStrings

#rename packages for convenience
dsm = DynamicStallModels
of = OpenFASTsr

# Get the path of the current file and change the working directory to that path
path = dirname(@__FILE__)
cd(path)

# Specify the airfoil and environmental variables, taken from pg 968
c = 1.5 # m, chord length
V = 60.0 #m/s, velocity 

# Read in the airfoil data 
Vertol = of.read_airfoilinput("../../data/airfoils/Vertol_23010-1.58.dat") #read in the airfoil data using OpenFASTsr
af = of.make_dsairfoil(Vertol) #make the airfoil into a DynamicStallModels airfoil
airfoils = Array{Airfoil, 1}(undef, 1) #make an array of the type Airfoil struct
airfoils[1] = af #put the airfoil into the array

# Make the Oye model struct
dsmodel = Oye(Indicial(), 1, airfoils) #makes the struct, says it will solve it indicially and tha there is 1 airfoil

# Create time, velocity, and angle of attack vectors
tvec = range(0, 0.05, 100) #time vector, these will be specific to the experimental data I am verifying against