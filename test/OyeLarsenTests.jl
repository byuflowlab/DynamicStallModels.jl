#=
OyeLarsenTests.jl
Jacob Child
Mar 25, 2023
Pseudocode: Use the Test package to run tests on the Oye/Larsen separation 
point functions and the Oye method overall
#Todo setup github to run this automatically on push
=#

#Packages to use
using Test, DynamicStallModels, DelimitedFiles, OpenFASTsr, Plots, CCBlade, FLOWMath 

dsm = DynamicStallModels
of = OpenFASTsr


#File to test if it is still good 
include("../testing/Oye/OyeComparer.jl")

#Come back to here 
path = dirname(@__FILE__)
cd(path)
#Setup to test the Oye/Larsen separation point functions (attachment degree f, and critical alphas (alphasep, alpha0?))

#Extract polar
Naca43618PolarTest = readdlm("../polars/LarsenNACA43618ClPolar.csv", ',') #read in the polar from Larsen
Naca43618PolarTest[:,1] = Naca43618PolarTest[:,1] .* pi/180 #convert to radians
Naca43618PolarTest = [Naca43618PolarTest zeros(size(Naca43618PolarTest[:,1]))]
#Extend polar
cr75 = .2 #chord/blade radius ratio at 75%
aoa, Cl, Cd = viterna(Naca43618PolarTest[:,1], Naca43618PolarTest[:,2], Naca43618PolarTest[:,3], cr75)
Cl[1] = Cl[2] - (Cl[3] - Cl[2]) #! fix this later? a cop out fix as the first value is really bad
polar = [aoa Cl Cd zeros(size(aoa))]
#Extract Fdyn data
FdynTest = readdlm("../testing/Oye/Outputs/LarsenFig4cFdynAttachment.csv", ',') #read in the Fdyn data from Larsen
#Quick for loop to make all Fdyn values close to 1 = 1 and values close to 0 = 0
for i in 1:length(FdynTest[:,2])
    if FdynTest[i,2] < .01
        FdynTest[i,2] = 0.0
    elseif FdynTest[i,2] > .99
        FdynTest[i,2] = 1.0
    end
    println(FdynTest[i,2])
end
FdynTestAk = Akima(FdynTest[:,1], FdynTest[:,2]) #make the akima interpolant
#Make airfoil struct
afTest = dsm.airfoil(polar; A = .07, sfun=dsm.LSP()) #make the airfoil struct
afTest = dsm.update_airfoil(afTest; alphasep=[afTest.alphasep[1], 32.0*pi/180] ) #update the airfoil with Larsen's separation point
airfoilsTest = Array{Airfoil, 1}(undef, 1) #make an array of the type Airfoil struct
airfoilsTest[1] = afTest #put the airfoil into the array
#Make the Oye struct and setup to solve
dsmodelTest = Oye(Indicial(), 1, airfoils,1,4) #makes the struct, says it will solve it indicially and that there is 1 airfoil

#Test the Oye/Larsen separation point functions (attachment degree f, and critical alphas (alphasep, alpha0?))

#= For reference:
fieldnames(typeof(afTest))
(:polar, :cl, :cd, :cm, :cn, :cc, :dcldalpha, :dcndalpha, :alpha0, :alphasep, :A, :b, :T, :sfun, :xcp, :eta, :zeta)
fieldnames(typeof(dsmodelTest))
(:detype, :n, :airfoils, :cflag, :version)
=#
#Test setup variables
alphaTestdeg = Vector(0:1.0:40.0)
TestTolerance = .01 #tolerance for the tests

@testset "Oye/Larsen Tests" begin 

    #Test set to test the attachment degree f
    @testset "f deg" begin
        @test dsm.separationpoint.(Ref(afTest), alphaTestdeg[3] .* pi/180) == FdynTestAk(alphaTestdeg[3]) #checking at about 2 deg, so should be fully attached
        @test dsm.separationpoint.(Ref(afTest), alphaTestdeg[11] .* pi/180) == FdynTestAk(alphaTestdeg[11]) #checking at about 10 deg, so linear region
        @test dsm.separationpoint.(Ref(afTest), alphaTestdeg[33] .* pi/180) == FdynTestAk(alphaTestdeg[33]) #checking at about 32 deg, so where deep stall begins
        @test dsm.separationpoint.(Ref(afTest), alphaTestdeg[41] .* pi/180) == FdynTestAk(alphaTestdeg[41]) #checking at about 40 deg, so in deep stall
        #print for debugging
        print("f af at 2 deg: ", dsm.separationpoint.(Ref(afTest), alphaTestdeg[2] .* pi/180), " fdynLarsen at 2 deg: ", FdynTestAk(alphaTestdeg[2] .* pi/180))
    end
    #Test set to test the critical alphas (alphasep, alpha0?)
    @testset "critical alphas" begin
        @test isapprox(afTest.alphasep, [-0.793985783199363, 0.5585053606381855 ], rtol = TestTolerance) #? should I do something about the low separation point at -45deg?
        @test isapprox(afTest.alpha0, -0.01982011981227349, rtol = TestTolerance) #! failing
        @test isapprox(afTest.dcldalpha, 6.634842198285429, rtol = TestTolerance) #! failing
        @test isapprox(afTest.dcndalpha, 6.634842198285429, rtol = TestTolerance) #! failing

    end
    #General Setup Tests to make sure the separation point function is Larsen's 
    @testset "Setup" begin 
        @test afTest.sfun == dsm.LSP()
        @test dsmodelTest.detype == Indicial()
        @test dsmodelTest.cflag == 1 #delay on lift 
        @test dsmodelTest.version == 4 #Larsen method
    end

end






#Plot the attachment degree f
#= How to do it with Vertol
dsm.separationpoint.(Ref(af), VertolPolar[:,1])
plot(VertolPolar[:,1] .*180/pi , dsm.separationpoint.(Ref(af), VertolPolar[:,1])
Note: to get plots that match my old f, I did not use LSP, I used adsp!
=#



#Plot the Larsen static polar and extrapolated polar to compare
plot(polar[:,1], polar[:,2], label = "Viterna Extrapolated Polar")
plot!(Naca43618PolarTest[:,1], Naca43618PolarTest[:,2], label="Larsen Static Polar")