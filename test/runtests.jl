using DynamicStallModels, DelimitedFiles, DifferentialEquations, OpenFASTTools, FLOWMath, Statistics
using Test

DE = DifferentialEquations
ds = DynamicStallModels
of = OpenFASTTools
@testset "DynamicStallModels" begin
    include("./utils_tests.jl") #Note: The out of domain warnings are typical... it's testing that it returns the correct value... outside of the domain. 
    include("./airfoil_tests.jl")
    include("./nomodel_tests.jl")
    include("./beddoesleishmanADG_tests.jl")
    include("./Oye_Test.jl")
    include("./Larsen_Separation_Point_Test.jl")
    include("Riso_Separation_Point_Tests.jl")
    include("riso_full_sep_tests.jl")
# include("./risotests.jl")
# include("./beddoesleishmantests.jl")
end

