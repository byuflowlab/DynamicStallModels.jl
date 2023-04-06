using DynamicStallModels, DelimitedFiles, DifferentialEquations, OpenFASTsr, FLOWMath, Statistics
using Test

DE = DifferentialEquations
ds = DynamicStallModels
of = OpenFASTsr

include("./utils_tests.jl")
include("./nomodel_tests.jl")
include("./beddoesleishmanADG_tests.jl")
# include("./risotests.jl")
# include("./beddoesleishmantests.jl")

