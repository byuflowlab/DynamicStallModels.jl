module DynamicStallModels

using FLOWMath, StaticArrays, Roots, DifferentialEquations, CurveFit

DE = DifferentialEquations

include("./utils.jl")
include("./types.jl")
include("./Riso/Riso.jl")
include("./BeddoesLeishman/BeddoesLeishman.jl")
include("./BeddoesLeishman/BeddoesLeishmanAeroDyn.jl")
include("./solve.jl")

end # module
