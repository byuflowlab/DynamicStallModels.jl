#=
Tests for all the utilities. 
=# 
using DynamicStallModels, DelimitedFiles, OpenFASTsr
using Test

DSM = DynamicStallModels
of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)

include("./testingutils.jl")

@testset "Utilities" begin
    @testset "Interpolation" begin
        ### Simple interpolation
        x = [0, 10]
        y = [0, 10]

        simpleinterp = Linear(x, y)

        #Test end points
        ygold = [0, 10]
        ytest = simpleinterp.(x)

        @test isapprox(ygold, ytest)
    
        #Test mid point
        ygold = 5.0
        ytest = simpleinterp(5.0)

        @test isapprox(ygold, ytest)

        #Test out of bounds
        xtest = [-1.0, 11]
        ytest = simpleinterp.(xtest)
        ygold = [0, 10]

        @test isapprox(ygold, ytest)


    end #End interpolation tests

    @testset "Utility Functions" begin
    end

    
end