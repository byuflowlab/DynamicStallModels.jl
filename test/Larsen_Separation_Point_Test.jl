using DynamicStallModels, DelimitedFiles, FLOWMath, DifferentialEquations, Test

dsm = DynamicStallModels

path = dirname(@__FILE__)
cd(path)

c = 0.55

M = 0.11
a = 343.0
Vrel = M*a 

file = "../polars/NACA_0015_Faber.csv"

polar = readdlm(file , ',')

dsmodel = Oye(Functional(), 1, 2, 4.0)

af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.LSP())

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af

tspan = (0, 2.0) #0:0.001:0.05

function Uvector(t)
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

parameters = [Uvector, 0.0, alpha, 0.0]

x_initial = [0.8]

prob = ODEProblem(airfoils, x_initial, tspan, parameters)

sol = DifferentialEquations.solve(prob, reltol=1e-8)

answer = parsesolution(dsmodel, af, sol, parameters)

@testset "Larsen Separation Point" begin
    @testset "Separation Point Values" begin
        #checks to see if the separation point is always between 0 and 1
        for i in 1:length(sol.t)
            Sep_Point = dsm.separationpoint(af, alpha(sol.t[i]))
            @test 0 <= Sep_Point <= 1
        end
    end
    @testset "Separation Angle of Attack" begin
        #checks to see if the separation angle of attacks are 32 degrees
        @test af.alphasep[2] == 32*pi/180
        @test af.alphasep[1] == -32*pi/180
    end
    @testset "Almost Fully Attached" begin
        #tests to see if close to the zero lift angle of attack the separation point is basically fully attached
        for i in af.alpha0:0.01:af.alpha0+(5*pi/180) 
            Sep_Point = dsm.separationpoint(af, i)
            @test isapprox(Sep_Point, 1.0, rtol = 0.05)
        end
    end
    @testset "Almost Fully Detatched" begin
        #tests to see if close to the fully separated angle of attack the separation point is close to 0.0
        for i in af.alphasep[2]-(5*pi/180):0.01:af.alphasep[2]
            Sep_Point = dsm.separationpoint(af, i)
            @test isapprox(Sep_Point, 0.0, atol = 0.05)
        end
    end
end
