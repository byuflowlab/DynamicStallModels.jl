module DynamicStallModels

using FLOWMath, StaticArrays, Roots, CurveFit


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

abstract type State end

include("./utils.jl")
include("./airfoils.jl")
include("./Oye.jl")
include("./Riso/Riso.jl")
include("./BeddoesLeishman/BeddoesLeishman.jl")
include("./BeddoesLeishman/BeddoesLeishmanAeroDyn.jl")
include("./BeddoesLeishman/BeddoesLeishmanADG.jl")
include("./solve.jl")

end # module
