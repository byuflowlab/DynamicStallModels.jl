#=
The no model model (for turning off Dynamic stall models but having a place holder). 
=#

struct NoModel <: DSModel
    detype::DEType
end

function NoModel()
    return NoModel(Discrete())
end

function numberofstates(dsmodel::NoModel)
    return 0
end

function initialize(dsmode::NoModel, airfoil::Airfoil, tvec, y, p)

    alpha = y[3]

    states = eltype(y)[]
    loads = [airfoil.cl(alpha), airfoil.cd(alpha), airfoil.cm(alpha)]

    return states, loads
end

function update_states!(model::NoModel, airfoil::Airfoil, x, xnew, y, p, dt)
end

function update_states(dsmodel::NoModel, airfoil::Airfoil, x, y, p, dt)
    return zeros(numberofstates(dsmodel))
end

function get_loads(dsmodel::NoModel, airfoil::Airfoil, states, y, p)
    alpha = y[3]
    cl = airfoil.cl(alpha)
    cd = airfoil.cd(alpha)
    cm = airfoil.cm(alpha)

    return cl, cd, cm
end

function get_loads!(dsmodel::NoModel, airfoil::Airfoil, states, loads, y, p)

    ### Retrieve the static loads. 
    alpha = y[3]
    cl = airfoil.cl(alpha) #Coefficient of lift
    cd = airfoil.cd(alpha) #Coefficient of drag
    cm = airfoil.cm(alpha) #Coefficient of moment

    #Store the loads in place. 
    loads[1] = cl
    loads[2] = cd
    loads[3] = cm
end