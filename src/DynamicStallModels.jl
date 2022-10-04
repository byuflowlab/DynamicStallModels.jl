module DynamicStallModels

using FLOWMath, StaticArrays, Roots, DifferentialEquations, CurveFit

DE = DifferentialEquations

export Functional, Iterative, Indicial

abstract type DEType end #Is the model designed to be solved in one go, or iteratively, updating p ever iteration. 

struct Functional <: DEType
end

struct Iterative <: DEType
end

struct Indicial <: DEType
end


abstract type DSModel end
# export DSModel

include("./utils.jl")
include("./types.jl")
include("./Riso/Riso.jl")
include("./BeddoesLeishman/BeddoesLeishman.jl")
include("./BeddoesLeishman/BeddoesLeishmanAeroDyn.jl")
include("./BeddoesLeishman/BeddoesLeishmanADG.jl")
include("./solve.jl")

end # module
