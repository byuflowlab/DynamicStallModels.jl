module DynamicStallModels

using FLOWMath, StaticArrays, Roots, CurveFit, SciMLBase, ForwardDiff


export Functional, Iterative, Indicial

abstract type DEType end #Is the model designed to be solved in one go, or iteratively, updating p ever iteration. 

struct Functional <: DEType #Todo: Rename to Continuous or something like that. 
end

struct Iterative <: DEType #Todo: Get rid of this. I don't think that I need it now that I figured out how to evaulate functions or iteratives within a functional. 
end

struct Indicial <: DEType #Todo: Rename to Discrete. 
end


abstract type DSModel end

const ny = 5


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
