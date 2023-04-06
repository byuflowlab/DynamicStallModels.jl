#=
Some interface functions for quick solving of the models for each airfoil. 

For the ode return, I can use methods on a struct on the airfoil. 

I can define a single method on a struct on the airfoil, that chains through to call whatever applicable update_states! or state_rates! function. Then for the vector of airfoils, you can define a method on that object, so 


Note: 
function (airfoil::Airfoil)(x, p, t)
end

function (airfoils::Array{Airfoil, 1})(x, p, t)
    return [airfoils[i].(x, p, t) for i in eachindex(airfoils)]
end

Or something like that. 
=#

export solve_indicial

function numberofstates_total(airfoils::Array{Airfoil,1})
    ns = 0

    for i in eachindex(airfoils)
        ns += numberofstates(airfoils[i].model)
    end
    return ns
end

function numberofparams_total(airfoils::Array{Airfoil,1})
    ns = 0

    for i in eachindex(airfoils)
        ns += numberofparams(airfoils[i].model)
    end
    return ns
end



function initialize(airfoil, tvec, y)  
    return initialize(airfoil.model, airfoil, tvec, y) 
end

function initialize_environment(U, udot, alpha, alphadot, n)
    return repeat([U, udot, alpha, alphadot], outer=n)
end

#=
    get_state_types(airfoils)

Get the different types of dynamic stall models used, and return a vector of their states so that they can be used to initialize an array of states. 
=#
function get_state_types(airfoils)
    n = length(airfoils)
    datatypes = Vector{DataType}(undef, n)
    for i = 1:n
        datatypes[i] = returnstatetype(airfoils[i].model)
    end
    return unique(datatypes)
end

function update_states(airfoil::Airfoil, oldstates, y, dt)
    newstates = zero(oldstates)
    return update_states!(airfoil.model, airfoil, oldstates, newstates, y, dt)
end

function update_states!(airfoil::Airfoil, oldstates, newstates, y, dt) #Note: Might need another function for state rate equations. 
    update_states!(airfoil.model, airfoil, oldstates, newstates, y, dt)
    return newstates
end

#=
Convinience function to access state rate equations
    - Todo: Add a check_dsm_type() function or capability to see if 
=#
function (airfoil::Airfoil)(x, p, t)
    if airfoil.model.detype==Indicial
        error("airfoil: This DEType cannot be used to return a state rate.")
    end

    state_rates!(airfoil.model, airfoil, dx, x, p, t)
end

function (airfoil::Airfoil)(dx, x, p, t)
end

function update_environment!(y, U, Udot, alpha, alphadot)
    y[1] = U
    y[2] = Udot
    y[3] = alpha
    y[4] = alphadot
end

function state_indices(dsmodel::DSModel, startidx)
    nsi = numberofstates(dsmodel)
    return startidx, startidx+nsi-1
end

function solve_indicial(airfoils::Array{Airfoil, 1}, tvec, Uvec, alphavec; verbose::Bool=false, Udotvec=zeros(length(Uvec)), alphadotvec=zeros(length(Uvec)))

    nt = length(tvec)
    n = length(airfoils)

    ### Initialize the states
    ns = numberofstates_total(airfoils)
    # np = numberofparams_total(airfoils)

    states = Array{eltype(Uvec), 2}(undef, nt, ns)
    y = initialize_environment(Uvec[1], Udotvec[1], alphavec[1], alphadotvec[1], n) #TODO: 
    stateidx = Vector{Int}(undef, n)
    loads = Array{eltype(Uvec), 2}(undef, nt, 3n) 

    tempx = 1
    # idxs = ns*(j-1)+1:j*ns 
    for i in eachindex(airfoils)
        stateidx[i] = tempx
        nsi1, nsi2 = state_indices(airfoils[i].model, stateidx[i])
        loadidx = 3*(i-1)+1:3i
        paramidx = 4*(i-1)+1:4*i #Todo: This needs to be updated to work with any number of input parameters -> I think I'm going to make the parameters be a constant length (the same for all DS models)
        ys = view(y, paramidx)
        states[1,nsi1:nsi2], loads[1,loadidx] = initialize(airfoils[i], tvec, ys) #Todo: make this function inplace. 
        tempx += numberofstates(airfoils[i].model)
    end

    # @show states[1,3]

    for i = 1:nt-1 
        ### Update environmental inputs
        t = tvec[i]
        dt = tvec[i+1]-tvec[i]

        

        if verbose
            @show t
        end

        for j = 1:n
            nsi1, nsi2 = state_indices(airfoils[j].model, stateidx[j])
            xsi = view(states, i, nsi1:nsi2)
            xs1 = view(states, i+1, nsi1:nsi2)
            ys = view(y, 4*(j-1)+1:4*j)

            idxs = 1:3
            loads_j = view(loads, i+1, idxs)

            update_environment!(ys, Uvec[i+1], Udotvec[i+1], alphavec[i+1], alphadotvec[i+1]) #TODO: Figure out how to make this work for varying stations. 
            
            update_states!(airfoils[j], xsi, xs1, ys, dt)

            getloads!(airfoils[j].model, airfoils[j], xs1, loads_j, ys) 
        end
    end

    return states, loads
end