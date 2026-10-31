using DynamicStallModels, DelimitedFiles, OpenFASTTools, FLOWMath, Statistics
using Test

#Note: DifferentialEquations is only needed by the commented-out test files below
#(Larsen_Separation_Point_Test.jl, risotests.jl, beddoesleishmantests.jl).
#Re-add `using DifferentialEquations` and `DE = DifferentialEquations` when those are enabled.
ds = DynamicStallModels
of = OpenFASTTools
@testset "DynamicStallModels" begin
    include("./utils_tests.jl") #Note: The out of domain warnings are typical... it's testing that it returns the correct value... outside of the domain. #Updated
    include("./airfoil_tests.jl") #Updated
    include("./nomodel_tests.jl") #Updated
    include("./beddoesleishmanADG_tests.jl") #Updated
    include("./Oye_Test.jl")
    include("./audit_regressions.jl")
    include("./repair_model_tests.jl")
    # include("./Larsen_Separation_Point_Test.jl") #Todo: Tests need to be cut down. 
# include("./risotests.jl")
# include("./beddoesleishmantests.jl")
end

