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

export solve_discrete, solve_indicial

function numberofstates_total(airfoils::AbstractVector{<:Airfoil})
    ns = 0

    for i in eachindex(airfoils)
        ns += numberofstates(airfoils[i].model) #todo: Same dynamic dispatch ... problem? I don't know if it's a problem. 
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



function initialize(airfoil, tvec, y, p)  #Todo: Why do I have tvec as one of the entries of this function? I didn't see something in the indicial models, so I wonder if I needed to initialize the continuous or discrete models. 
    return initialize(airfoil.model, airfoil, tvec, y, p) 
end

function initialize_environment(U, udot, alpha, alphadot, n)
    return repeat([U, udot, alpha, alphadot], outer=n)
end

# The unused private get_state_types helper was removed: state storage depends
# on the supplied numerical inputs, not just on a model's type.

function update_states(airfoil::Airfoil, oldstates, y, p, dt)
    newstates = similar(oldstates,promote_type(eltype(oldstates),eltype(airfoil.polar)))
    return update_states!(airfoil.model, airfoil, oldstates, newstates, y, p, dt)
end

function update_states!(airfoil::Airfoil, oldstates, newstates, y, p, dt) #Note: Might need another function for state rate equations. 
    update_states!(airfoil.model, airfoil, oldstates, newstates, y, p, dt)
    return newstates
end

function state_rates!(airfoil, state_rates, states, y, p, t_aspect)
    # @show state_rates
    # @show states
    state_rates!(airfoil.model, airfoil, state_rates, states, y, p, t_aspect)
    # @show state_rates
    # @show states
    # println("")
end

function get_loads(dsmodel::DSModel, states, y, p, airfoil)
    loads = zeros(3)
    get_loads!(dsmodel, airfoil, states, loads, y, p)
    return loads
end

#=
Convinience function to access state rate equations
    - Todo: Add a check_dsm_type() function or capability to see if 
=#
function (airfoil::Airfoil)(state_in, y, p, t_aspect)
    state_out = similar(state_in,promote_type(eltype(state_in),eltype(airfoil.polar)),numberofstates(airfoil.model))
    return airfoil(state_out, state_in, y, p, t_aspect)
end

# Same as above, but for calculations
"""
    (airfoil::Airfoil)(state_out, state_in, y, p, t_aspect)

A method on the Airfoil struct to either update the states, or return the state rates of the dynamic stall model (respective to if the model is discrete, or continuous). 

### Inputs
- state_in::Vector{TF} - A vector of the old states. 
- state_out::Vector{TF} - A vector of the new states if the model is discrete, and a vector of the new state rates if the model is a state space model (continuous). 
- y::Vector{TF} - A vector of the environmental inputs. 
- t_aspect::TF - Either the time step if the model is discrete, or the current time value if the model is a state space model. 
"""
function (airfoil::Airfoil)(state_out, state_in, y, p, t_aspect)
    if isa(airfoil.model.detype, Discrete)
        return update_states!(airfoil, state_in, state_out, y, p, t_aspect)
    elseif isa(airfoil.model.detype, Continuous)
        return state_rates!(airfoil, state_out, state_in, y, p, t_aspect)
    end
end

function (airfoils::AbstractVector{<:Airfoil})(state_idxs, state_in, y, p, t_aspect)
    ns = numberofstates_total(airfoils)
    state_out = similar(state_in,promote_type(eltype(state_in),map(af->eltype(af.polar),airfoils)...),ns)
    airfoils(state_out, state_in, state_idxs, y, p, t_aspect)
    return state_out
end


function (airfoils::AbstractVector{<:Airfoil})(state_out, state_in, state_idxs, y, p, t_aspect)
    n=length(airfoils)
    length(state_idxs)==n || throw(DimensionMismatch("Provide a state offset per airfoil."))
    length(state_in)==length(state_out)==numberofstates_total(airfoils) || throw(DimensionMismatch("State sizes must match the airfoils."))
    length(y)==4n && length(p)==2n || throw(DimensionMismatch("Provide four environment and two geometry entries per airfoil."))
    
    for i in eachindex(airfoils)
        nsi1, nsi2 = state_indices(airfoils[i].model, state_idxs[i])
        xsi = view(state_in, nsi1:nsi2)
        xs1 = view(state_out, nsi1:nsi2)
        ys = view(y, 4*(i-1)+1:4*i)
        ps = view(p, 2*(i-1)+1:2*i)

        airfoils[i](xs1, xsi, ys, ps, t_aspect)
    end
    return state_out
end

function update_environment!(y, U, Udot, alpha, alphadot)
    y[1] = U
    y[2] = Udot
    y[3] = alpha
    y[4] = alphadot
end

function state_indices(dsmodel::DSModel, startidx) #todo: -> should this function be dispatching on the type of dsmodel? Like have a state_indices function for each each dsmodel? 
    nsi = numberofstates(dsmodel) #todo: nsi is assigned as Union{Nothing, Int64} -> numberofstates() probably needs to require exact outputs? 
    return startidx, startidx+nsi-1 #todo: dynamic dispatch to DSM.numberofstates()
end

"""
    solve_discrete(airfoils, t, U, alpha; cvec=ones(length(airfoils)),
                   xcpvec=fill(0.2, length(airfoils)),
                   Udotvec=zeros(length(U)), alphadotvec=zeros(length(U)), verbose=false)

Advance Discrete() models over a finite, strictly increasing time grid (s).
The sampled environment is shared: U in m/s, alpha in radians, Udot in m/s²,
alphadot in rad/s. Chords are in metres and xcp is nondimensional. Each step
holds its right-end input constant. Return time-by-state states and loads
grouped as Cl, Cd, Cm per airfoil. The first row is the model initialization.
solve_indicial remains an exact compatibility alias.
"""
function solve_discrete(airfoils::AbstractVector{<:Airfoil}, tvec, Uvec, alphavec; verbose::Bool=false, cvec=ones(length(airfoils)), xcpvec=0.2.*ones(length(airfoils)), Udotvec=zeros(length(Uvec)), alphadotvec=zeros(length(Uvec)))

    nt = length(tvec)
    n = length(airfoils)
    n > 0 && nt > 0 || throw(ArgumentError("Provide at least one airfoil and time sample."))
    all(isfinite, tvec) && all(diff(tvec) .> 0) || throw(ArgumentError("Time must be finite and strictly increasing."))
    all(length(v) == nt for v in (Uvec, alphavec, Udotvec, alphadotvec)) || throw(DimensionMismatch("Environment histories must match the time grid."))
    length(cvec) == length(xcpvec) == n || throw(DimensionMismatch("Provide one chord and xcp per airfoil."))
    all(af.model.detype isa Discrete for af in airfoils) || throw(ArgumentError("solve_discrete requires Discrete() models."))
    all(v->all(isfinite,v),(Uvec,alphavec,Udotvec,alphadotvec,cvec,xcpvec)) || throw(ArgumentError("Environment and geometry must be finite."))
    all(Uvec .>= 0) || throw(ArgumentError("Relative-speed magnitudes must be nonnegative."))
    all(cvec .> 0) || throw(ArgumentError("Chord lengths must be positive."))

    ### Initialize the states
    ns = numberofstates_total(airfoils)
    # np = numberofparams_total(airfoils)
    inittype = find_inittype(cvec[1], Uvec[1], alphavec[1])

    states = Array{inittype, 2}(undef, nt, ns)
    y = initialize_environment(Uvec[1], Udotvec[1], alphavec[1], alphadotvec[1], n) 
    # p = cvec

    stateidx = Vector{Int}(undef, n)
    p = Vector{inittype}(undef, 2n)
    loads = Array{inittype, 2}(undef, nt, 3n) 

    tempx = 1
    # idxs = ns*(j-1)+1:j*ns 
    for i in eachindex(airfoils)
        stateidx[i] = tempx
        nsi1, nsi2 = state_indices(airfoils[i].model, stateidx[i])
        loadidx = 3*(i-1)+1:3i
        envidx = 4*(i-1)+1:4*i #todo: This needs to be updated to work with any number of input parameters -> I think I'm going to make the parameters be a constant length (the same for all DS models)
        paramidx = 2*(i-1)+1:2*i
        ys = view(y, envidx)

        p[paramidx] = [cvec[i], xcpvec[i]] 
        ps = view(p, paramidx)

        states[1,nsi1:nsi2], loads[1,loadidx] = initialize(airfoils[i], tvec, ys, ps) #todo: make this function inplace. 
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
            ps = view(p, 2*(j-1)+1:2*j)
         
            
            idxs = 3*(j-1)+1:3j
            loads_j = view(loads, i+1, idxs)

            update_environment!(ys, Uvec[i+1], Udotvec[i+1], alphavec[i+1], alphadotvec[i+1]) #TODO: Figure out how to make this work for varying stations. -> I don't know what I meant by that. 
            
            airfoils[j](xs1, xsi, ys, ps, dt)

            get_loads!(airfoils[j].model, airfoils[j], xs1, loads_j, ys, ps) 
        end
    end

    return states, loads
end

const solve_indicial = solve_discrete

export ODEProblem

"""
    ODEProblem(airfoils, x0, tspan, y; cvec=ones(length(airfoils)),
               xcpvec=fill(0.2, length(airfoils)))

Construct an in-place SciML problem for Continuous() models. For each airfoil,
y contains U, Udot, alpha, alphadot (numbers or functions of time). Geometry
is supplied separately as one chord and xcp per airfoil.
"""
function SciMLBase.ODEProblem(airfoils::AbstractVector{<:Airfoil}, x0, tspan, y; cvec=ones(length(airfoils)), xcpvec=fill(0.2, length(airfoils)))
    #Todo: account for x, y, p, t (p being paramters). -> GXBeam probably has some examples. 

    n = length(airfoils)
    length(y) == 4n || throw(DimensionMismatch("Provide four environment entries per airfoil."))
    length(cvec) == length(xcpvec) == n || throw(DimensionMismatch("Provide one chord and xcp per airfoil."))
    length(x0) == numberofstates_total(airfoils) || throw(DimensionMismatch("Initial state count does not match the airfoils."))
    all(af.model.detype isa Continuous for af in airfoils) || throw(ArgumentError("ODEProblem requires Continuous() models."))
    geometry = collect(Iterators.flatten(zip(cvec, xcpvec)))
    stateidx = Vector{Int}(undef, n)

    tempx = 1
    
    for i in eachindex(airfoils)
        stateidx[i] = tempx
        tempx += numberofstates(airfoils[i].model)
    end

    fun(dx, x, p, t) = airfoils(dx, x, stateidx, p, geometry, t)

    # @show typeof(tspan)

    return SciMLBase.ODEProblem{true, true}(fun, x0, tspan, y)
end

function evaluate_environment(y, t)

    U = isa(y[1], Function) ? y[1](t) : y[1]
    Udot = isa(y[2], Function) ? y[2](t) : y[2]
    alpha = isa(y[3], Function) ? y[3](t) : y[3]
    alphadot = isa(y[4], Function) ? y[4](t) : y[4]

    return U, Udot, alpha, alphadot
end
