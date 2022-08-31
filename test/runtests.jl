using DynamicStallModels, DelimitedFiles, DifferentialEquations, OpenFASTsr, FLOWMath, Statistics
using Test

DE = DifferentialEquations
ds = DynamicStallModels
of = OpenFASTsr

include("./risotests.jl")
include("./beddoesleishmantests.jl")