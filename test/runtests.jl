using DynamicStallModels, DelimitedFiles, DifferentialEquations, OpenFASTTools, FLOWMath, Statistics
using Test

DE = DifferentialEquations
ds = DynamicStallModels
of = OpenFASTTools
@testset "DynamicStallModels" begin
    include("./utils_tests.jl") #Note: The out of domain warnings are typical... it's testing that it returns the correct value... outside of the domain. #Updated
    include("./airfoil_tests.jl") #Updated
    include("./nomodel_tests.jl") #Updated
    include("./beddoesleishmanADG_tests.jl") #Updated
    # include("./Oye_Test.jl") #Todo: Tests need to be trimmed down. -> There are failing tests. 
    # include("./Larsen_Separation_Point_Test.jl") #Todo: Tests need to be cut down. 
# include("./risotests.jl")
# include("./beddoesleishmantests.jl")
end

