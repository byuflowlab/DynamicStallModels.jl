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

#Come back to here 
path = dirname(@__FILE__)
cd(path)
#Setup to test the Oye/Larsen separation point functions (attachment degree f, and critical alphas (alphasep, alpha0?))

#Extract polar
Naca43618PolarTest = readdlm("../polars/LarsenFig4bClStaticPlot.csv", ',') #read in the polar from Larsen
Naca43618PolarTest[:,1] = Naca43618PolarTest[:,1] .* pi/180 #convert to radians
Naca43618PolarTest = [Naca43618PolarTest zeros(size(Naca43618PolarTest[:,1]))]
#Extend polar
#TODO before extending add a point at [0,0]
cr75 = .2 #chord/blade radius ratio at 75%
aoa, Cl, Cd = viterna(Naca43618PolarTest[:,1], Naca43618PolarTest[:,2], Naca43618PolarTest[:,3], cr75)
#Cl[1] = 0.0 #Cl[2] - (Cl[3] - Cl[2]) #! fix this later? a cop out fix as the first value is really bad
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
    #println(FdynTest[i,2])
end


#Make airfoil struct
#! somehow the airfoil function is being called twice?
afTest = dsm.airfoil(polar; A = .07, sfun=dsm.LSP()) #make the airfoil struct
afTest = dsm.update_airfoil(afTest; alphasep=[afTest.alphasep[1], 30.0*pi/180], dcldalpha = 2* pi, dcndalpha = 2 * pi ) #update the airfoil with Larsen's separation point

#! I deleted FdynTestAk as we should compare just to explicit data
#for the sake of testing ->
delta = zeros(length(FdynTest[:,1]) - 1)
for i in 1:(length(FdynTest[:,1])-2)
    delta[i] = FdynTest[i+1,1] - FdynTest[i,1]
end
avg = sum(delta)/length(delta)
FdynExtraAlpha = Vector(afTest.alpha0-10*avg:avg:FdynTest[1,1])
FdynExtraF = zeros(length(FdynExtraAlpha))
FFull = [[FdynExtraAlpha FdynExtraF]; FdynTest]
#Find where FFull starts going down
FirstOne = findfirst(x -> x == 1.0, FFull[:,2])
println("FFull starts going back down at: ", FFull[findlast(x -> x == 1.0, FFull[:,2]),1], " deg")

#Plot and compare the F from dsm to Fdyn from Larsen
plot(FFull[:,1], FFull[:,2], label = "Larsen", title = "Fdyn Comparison", xlabel = "alpha (deg)", ylabel = "Fdyn")
plot!(FFull[:,1], dsm.separationpoint.(Ref(afTest), FFull[:,1] .* pi/180,avg), label = "DSM")


#=
#! the following code is just to be ran once and will not work again
cnTest = zeros(length(FFull[:,1]))
cn_sepTest = zeros(length(FFull[:,1]))
cn_invTest = zeros(length(FFull[:,1]))
cn_fsTest = zeros(length(FFull[:,1]))
fstTest = zeros(length(FFull[:,1]))
for i in 1:20
    cnTest[i], cn_sepTest[i], cn_invTest[i], cn_fsTest[i], fstTest[i] = dsm.separationpoint.(Ref(afTest), FFull[i,1] .* pi/180)
end
#plot all the data close up 
plot(FFull[1:20,1], cnTest[1:20], label = "Cl static", title = "Cl Comparison", xlabel = "alpha (deg)", ylabel = "Cl")
plot!(FFull[1:20,1], cn_sepTest[1:20], label = "Cl sep")
plot!(FFull[1:20,1], cn_invTest[1:20], label = "Cl inv")
plot!(FFull[1:20,1], cn_fsTest[1:20], label = "Cl fs", legend = :bottomright)
plot!(FFull[1:20,1], FFull[1:20,2], label = "F Larsen")
plot!(FFull[1:20,1], dsm.separationpoint.(Ref(afTest),FFull[1:20,1] .*pi/180), label = "F DSM")
plot!(ylabel = "Cl/Fdyn")
=#

#make my AlphaTest vector
#AlphaTest = Vector(-20:0.1:40.0)


#= #commenting all out so I can figure out the if statement fdyn follower
airfoilsTest = Array{Airfoil, 1}(undef, 1) #make an array of the type Airfoil struct
airfoilsTest[1] = afTest #put the airfoil into the array
#Make the Oye struct and setup to solve
dsmodelTest = Oye(Indicial(), 1, airfoilsTest,1,4) #makes the struct, says it will solve it indicially and that there is 1 airfoil

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

#Plot the Fdyn comparision
plot(alphaTestdeg, dsm.separationpoint.(Ref(afTest), alphaTestdeg.*pi/180), label = "DSM")
plot!(alphaTestdeg, FdynTestAk(alphaTestdeg), label = "Larsen")


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
    @testset "critical alphas" begin #! where did these values come from?
        @test isapprox(afTest.alphasep, [-0.793985783199363, 0.5585053606381855 ], rtol = TestTolerance) #? should I do something about the low separation point at -45deg?
        @test isapprox(afTest.alpha0, -0.01982011981227349, rtol = TestTolerance) #! failing
        @test isapprox(afTest.dcldalpha, 2*pi, rtol = TestTolerance) #! failing
        @test isapprox(afTest.dcndalpha, 2*pi, rtol = TestTolerance) #! failing

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
#=
plot(polar[:,1], polar[:,2], label = "Viterna Extrapolated Polar")
plot!(Naca43618PolarTest[:,1], Naca43618PolarTest[:,2], label="Larsen Static Polar")
=#


=#