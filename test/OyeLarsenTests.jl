#=
OyeLarsenTests.jl
Jacob Child
Mar 25, 2023
Pseudocode: Use the Test package to run tests on the Oye/Larsen separation 
point functions and the Oye method overall
#Todo setup github to run this automatically on push
=#

#Packages to use
using Test, DynamicStallModels, DelimitedFiles, OpenFASTsr, Plots 

DSM = DynamicStallModels
of = OpenFASTsr

#Setup to test the Oye/Larsen separation point functions (attachment degree f, and critical alphas (alphasep, alpha0?))
Naca43618Polar = readdlm("../polars/LarsenNACA43618ClPolar.csv", ',') #read in the polar from Larsen
af = dsm.airfoil(Naca43618Polar; A = .07, sfun=dsm.LSP()) #make the airfoil struct
af = dsm.update_airfoil(af; alphasep=[af.alphasep[1], 32.0*pi/180] ) #update the airfoil with Larsen's separation point
airfoils = Array{Airfoil, 1}(undef, 1) #make an array of the type Airfoil struct
airfoils[1] = af #put the airfoil into the array
dsmodel = Oye(Indicial(), 1, airfoils,1,4) #makes the struct, says it will solve it indicially and that there is 1 airfoil


#Plot the attachment degree f
#= How to do it with Vertol
dsm.separationpoint.(Ref(af), VertolPolar[:1])
plot(VertolPolar[:,1] .*180/pi , dsm.separationpoint.(Ref(af), VertolPolar[:,1])
Note: to get plots that match my old f, I did not use LSP, I used adsp!

=#