module DynamicStallModels

using FLOWMath, StaticArrays, Roots, CurveFit, SciMLBase, ForwardDiff

#Todo. Removed Roots for compatibility with Julia 1.11. Affects Larsen and Riso models. 


export Continuous, Discrete

abstract type DEType end #Is the model designed to be solved in one go, or iteratively, updating p ever iteration. 

struct Continuous <: DEType 
end

struct Discrete <: DEType
end


abstract type DSModel end


include("./utils.jl")
include("./airfoils.jl")
include("./NoModel.jl")
include("./Oye.jl")
include("./Riso/Riso.jl")
include("./BeddoesLeishman/BeddoesLeishman.jl")
include("./BeddoesLeishman/BeddoesLeishmanAeroDyn.jl")
include("./BeddoesLeishman/BeddoesLeishmanADG.jl")
include("./solve.jl")

end # module
