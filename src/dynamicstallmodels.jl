module dynamicstallmodels

using FLOWMath, StaticArrays, Roots, DifferentialEquations, CurveFit

DE = DifferentialEquations

include("./utils.jl")
include("./types.jl")
include("./Riso/Riso.jl")

end # module
