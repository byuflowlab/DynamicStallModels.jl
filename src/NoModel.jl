#=
The no model model (for turning off Dynamic stall models but having a place holder). 
=#

struct NoModel <: DSModel
    detype::DEType
end

function NoModel()
    return NoModel(Indicial())
end

function numberofstates(dsmodel::NoModel)
    return 0
end

function initialize(dsmode::NoModel, airfoil::Airfoil, tvec, y)

    states = eltype(y)[]
    loads = zeros(3)

    return states, loads
end

function update_states!(model::NoModel, airfoil::Airfoil, x, xnew, y, dt)
end

function update_states(dsmodel::NoModel, airfoil::Airfoil, x, y, dt)
    return zeros(numberofstates(dsmodel))
end

function get_loads!(dsmodel::NoModel, airfoil::Airfoil, states, loads, y)
    #Todo: This should just grab a load from the airfoil polar. 
    # println("Got here")
    loads[:] .= 0.0
end