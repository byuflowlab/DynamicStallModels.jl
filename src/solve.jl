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

function numberofstates_total(airfoils::AbstractVector{<:Airfoil})
    ns = 0

    for i in eachindex(airfoils)
        ns += numberofstates(airfoils[i].model)
    end
    return ns
end

function numberofparams_total(airfoils::AbstractVector{<:Airfoil})
    ns = 0

    for i in eachindex(airfoils)
        ns += numberofparams(airfoils[i].model)
    end
    return ns
end



function initialize(airfoil, tvec, y)  #Todo: Why do I have tvec as one of the entries of this function? I didn't see something in the indicial models, so I wonder if I needed to initialize the functionals or iterative models. 
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

function get_loads(dsmodel::DSModel, states, y, airfoil)
    loads = zeros(3)
    get_loads!(dsmodel, airfoil, states, loads, y)
end

#=
Convinience function to access state rate equations
    - Todo: Add a check_dsm_type() function or capability to see if 
=#
function (airfoil::Airfoil)(state_in, y, t_aspect)
    state_out = zeros(numberofstates(airfoil.model))
    return airfoil(state_out, state_in, y, t_aspect)
end

# Same as above, but for calculations
"""
    (airfoil::Airfoil)(state_in, state_out, y, t_aspect)

A method on the Airfoil struct to either update the states, or the state rates of the dynamic stall model (respective to if the model is indicial, or functional/iterative). 

### Inputs
- state_in::Vector{TF} - A vector of the old states. 
- state_out::Vector{TF} - A vector of the new states if the model is indicial, and a vector of the new state rates if the model is a state space model (functional/iterative). 
- y::Vector{TF} - A vector of the environmental inputs. 
- t_aspect::TF - Either the time step if the model is indicial, or the current time value if the model is a state space model. 
"""
function (airfoil::Airfoil)(state_out, state_in, y, t_aspect)
    if isa(airfoil.model.detype, Indicial)
        return update_states!(airfoil, state_out, state_in, y, t_aspect)
    elseif isa(airfoil.model.detype, Functional)
        return state_rates!(airfoil, state_out, state_in, y, t_aspect)
    end
end

function (airfoils::AbstractVector{<:Airfoil})(state_idxs, state_in, y, t_aspect)
    ns = numberofstates_total(airfoils)
    state_out = zeros(ns)

    for i in eachindex(airfoils)
        nsi1, nsi2 = state_indices(airfoils[i].model, state_idxs[i])
        xsi = view(state_out, trunc(Int,nsi1+1):trunc(Int,nsi2+1))
        xs1 = view(state_in, trunc(Int,nsi1+1):trunc(Int,nsi2+1))
        ys = view(y, 4*(i-1)+1:4*i)


        airfoils[i](xsi, xs1, ys, t_aspect)

        # if any(item -> isnan(item), xs1) #Todo: I don't know if having something like this would be more robust, or just a nusiance. 
        #     @show xsi
        #     @show xs1
        #     @show ys
        #     error("Found NaN while updating states")
        # end
    end
end


function (airfoils::AbstractVector{<:Airfoil})(state_in, state_out, state_idxs, y, t_aspect)

    for i in eachindex(airfoils)
        nsi1, nsi2 = state_indices(airfoils[i].model, state_idxs[i])
        xsi = view(state_in, nsi1:nsi2)
        xs1 = view(state_out, nsi1:nsi2)
        ys = view(y, 4*(i-1)+1:4*i)

        # @show ys

        airfoils[i](xsi, xs1, ys, t_aspect)
    end
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

            # idxs = 1:3 #Todo: What the heck is going on here? 
            idxs = 3*(j-1)+1:3j
            loads_j = view(loads, i+1, idxs)

            update_environment!(ys, Uvec[i+1], Udotvec[i+1], alphavec[i+1], alphadotvec[i+1]) #TODO: Figure out how to make this work for varying stations. 
            
            airfoils[j](xsi, xs1, ys, dt) #Todo: I had xsi, xs1... 

            get_loads!(airfoils[j].model, airfoils[j], xs1, loads_j, ys) 
        end
    end

    return states, loads
end



