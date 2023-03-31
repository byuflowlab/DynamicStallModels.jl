#=
Some interface functions for quick solving of the models for each airfoil. 
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



function initialize(airfoil, tvec, U, alpha)  
    return initialize(airfoil.model, airfoil, tvec, U, alpha) 
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

function state_indices(dsmodel::DSModel, startidx)
    nsi = numberofstates(dsmodel)
    return startidx, startidx+nsi-1
end

function solve_indicial(airfoils::Array{Airfoil, 1}, tvec, Uvec, alphavec; verbose::Bool=false)

    nt = length(tvec)
    n = length(airfoils)

    ### Initialize the states
    ns = numberofstates_total(airfoils)
    np = numberofparams_total(airfoils)

    states = Array{eltype(Uvec), 2}(undef, nt, ns)
    y = Vector{eltype(Uvec)}(undef, np)
    stateidx = Vector{Int}(undef, n)
    loads = Array{eltype(Uvec), 2}(undef, nt, 3n) 

    tempx = 1
    # idxs = ns*(j-1)+1:j*ns 
    for i in eachindex(airfoils)
        stateidx[i] = tempx
        nsi1, nsi2 = state_indices(airfoils[i].model, stateidx[i])
        loadidx = 3*(i-1)+1:3i
        paramidx = 2*(i-1)+1:2i #Todo: This needs to be updated to work with any number of input parameters -> I think I'm going to make the parameters be a constant length (the same for all DS models)
        states[1,nsi1:nsi2], loads[1,loadidx], y[paramidx] = initialize(airfoils[i], tvec, Uvec[1], alphavec[1]) #Todo: make this function inplace. 
        tempx += numberofstates(airfoils[i].model)
    end

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
            ys = view(y, 2*(j-1)+1:2j)

            idxs = 1:3
            loads_j = view(loads, i+1, idxs)

            update_environment!(airfoils[j].model, ys, Uvec[i+1], alphavec[i+1]) #Todo: This probably can be a universal function. 
            
            update_states!(airfoils[j], xsi, xs1, ys, dt)

            getloads!(airfoils[j].model, airfoils[j], xs1, loads_j, ys) 
        end
    end

    return states, loads
end