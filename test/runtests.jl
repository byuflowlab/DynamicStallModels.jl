using DynamicStallModels, DelimitedFiles, DifferentialEquations, OpenFASTsr, FLOWMath, Statistics
using Test

DE = DifferentialEquations
ds = DynamicStallModels
of = OpenFASTsr
@testset "DynamicStallModels" begin
    include("./utils_tests.jl") #Note: The out of domain warnings are typical... it's testing that it returns the correct value... outside of the domain. 
    include("./airfoil_tests.jl")
    include("./nomodel_tests.jl")
    include("./beddoesleishmanADG_tests.jl")
    # include("./Oye_Test.jl") #Todo: Tests need to be cut down. -> There are failing tests. 
    # include("./Larsen_Separation_Point_Test.jl") #Todo: Tests need to be cut down. 
# include("./risotests.jl")
# include("./beddoesleishmantests.jl")
end

