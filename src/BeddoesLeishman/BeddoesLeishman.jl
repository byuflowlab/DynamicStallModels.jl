#=
The state space form of the Beddoes-Leishman model. Including several different implementations. 


Adam Cardoza 8/24/22
=#

export BeddoesLeishman

"""
    BeddoesLeishman(detype::DEType, n::Int, airfoils::Array{Airfoil, 1}, version::Int)

The Beddoes-Leishman model struct. It stores airfoil data for every section to be simulated. It can be used as a method to return updated states or state rates depending on it's DEType. 

### Inputs
- detype - The type of model it is, Functional(), Iterative(), or Indicial().
- n - The number of sections to be simulated. 
- airfoils - A vector of Airfoil structs, one corresponding to each section to be simulated. 
- version - Which version of the indicial implementation. 1) original. 2) AeroDyn original. 3) AeroDyn Gonzalez. 4) AeroDyn Minema
"""
struct BeddoesLeishman{TI} <: DSModel
    detype::DEType 
    n::TI #Number of airfoils simulated
    airfoils::Array{Airfoil,1} #Todo: This is inconvinient typing.
    version::TI #Which version of the indicial implementation.  
end

function (model::BeddoesLeishman)(x, p, t, dt) 
    if isa(model.detype, Functional)
        @warn("Functional implementation not yet prepared.")
    elseif isa(model.detype, Indicial)
        if model.version==1
            @warn("Original indicial Beddoe-Leishman not prepared for use yet.")
        elseif model.version==2
            ns = numberofstates(model)
            newstates = Array{eltype(p), 1}(undef, ns)
            for i = 1:model.n
                ps = view(p, 18*(i-1)+1:18*i)
                #[c, a, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, xcp]
                c, a, dcndalpha, alpha0, _, _, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, _, U, aoa = ps #Inputs  
                
                xs = view(x, 22*(i-1)+1:22*i)

                idx = 22*(i-1)+1:22*i
                newstates[idx] = update_states_ADO(model, xs, c, a, U, dt, aoa, dcndalpha, alpha0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, i)
            end
            return newstates

        elseif model.version==3
            # @warn("AeroDyn Beddoe-Leishman with Gonzalez's modifications not prepared for use yet.")

            nst = numberofstates_total(model)
            ns = numberofstates(model)
            np = numberofparams(model)
            newstates = Array{eltype(p), 1}(undef, nst)
            for i = 1:model.n
                ps = view(p, np*(i-1)+1:np*i)
    
                c, a, dcndalpha, alpha0, _, _, A1, A2, A5, b1, b2, b5, Tf0, Tv0, Tp, Tvl, Tsh, Cn1, _, zeta, U, aoa = ps #Inputs  
                
                xs = view(x, ns*(i-1)+1:ns*i)

                idx = ns*(i-1)+1:ns*i

                newstates[idx] = update_states_ADG(model, xs, c, a, U, dt, aoa, dcndalpha, alpha0, A1, A2, A5, b1, b2, b5, Tf0, Tv0, Tp, Tvl, Tsh, Cn1, zeta, i)
            end
            return newstates

        elseif model.version==4
            @warn("AeroDyn Beddoe-Leishman with Minema's modifications not prepared for use yet.")
        end
    end
end

function getloads(dsmodel::BeddoesLeishman, states, p, airfoil)
    if isa(dsmodel.detype, Functional)
        @warn("Functional implementation not yet prepared.")
    elseif isa(dsmodel.detype, Indicial)
        if dsmodel.version==1
            @warn("Original indicial Beddoe-Leishman not prepared for use yet.")
        elseif dsmodel.version==2
            return getloads_BLA(dsmodel, states, p, airfoil)
        elseif dsmodel.version==3
            # @warn("AeroDyn Beddoe-Leishman with Gonzalez's modifications not prepared for use yet.")
            return getloads_BLAG(dsmodel, states, p, airfoil)
        elseif dsmodel.version==4
            @warn("AeroDyn Beddoe-Leishman with Minema's modifications not prepared for use yet.")
        end
    end
end

function numberofstates(dsmodel::BeddoesLeishman) #TODO: This probably need to be augmented to check if the model is a functional, an iterative, or an indicial. 
    if dsmodel.version==1 #Todo: This function should be changed to reflect the number of states for a single 2D airfoil, not the total number of states.... or I should have two different functions. 
        @warn("The orginal Beddoes-Leishman model is not yet prepared.")
        return 0
    elseif dsmodel.version==2
        return 22
    elseif dsmodel.version==3
        # @warn("The Gozalez Beddoes-Leishman model is not yet prepared.")
        return 32
    elseif dsmodel.version==4
        @warn("The Minema Beddoes-Leishman model is not yet prepared.")
        return 22
    end
end

function numberofloads(dsmodel::BeddoesLeishman) #TODO: This probably need to be augmented to check if the model is a functional, an iterative, or an indicial.
    if dsmodel.version==1
        @warn("The orginal Beddoes-Leishman model is not yet prepared.")
        return 0
    elseif dsmodel.version==2
        return 5
    elseif dsmodel.version==3
        # @warn("The Gonzalez Beddoes-Leishman model is not yet prepared.")
        return 5
    elseif dsmodel.version==4
        @warn("The Minema Beddoes-Leishman model is not yet prepared.")
        return 5
    end
end

function numberofparams(dsmodel::BeddoesLeishman) #TODO: This probably need to be augmented to check if the model is a functional, an iterative, or an indicial.

    if dsmodel.version==1
        @warn("The orginal Beddoes-Leishman model is not yet prepared.")
        return 0
    elseif dsmodel.version==2
        return 22
    elseif dsmodel.version==3
        # @warn("The Gonzalez Beddoes-Leishman model is not yet prepared.")
        return 22
    elseif dsmodel.version==4
        @warn("The Minema Beddoes-Leishman model is not yet prepared.")
        return 22
    end
end

function initialize(dsmodel::BeddoesLeishman, Uvec, aoavec, tvec, airfoil::Airfoil, c, a)
    if isa(dsmodel.detype, Functional)
        @warn("Beddoes Leishman Functional implementation isn't prepared yet. - initialize()")
    elseif isa(dsmodel.detype, Iterative)
        @warn("Beddoes Leishman Iterative implementation isn't prepared yet. - initialize()")
    else #Model is indicial
        if dsmodel.version==2 # Original AeroDyn Implementation
            return initialize_ADO(Uvec, aoavec, tvec, airfoil, c, a)
        elseif dsmodel.version==3 # AeroDyn with Gonzalez modifications. 
            return initialize_ADG(Uvec, aoavec, tvec, airfoil, c, a)
        end
    end
end

function update_environment!(dsmodel::BeddoesLeishman, p, U, aoa)
    if isa(dsmodel.detype, Functional)
        @warn("Beddoes Leishman Functional implementation isn't prepared yet. - initialize()")
    elseif isa(dsmodel.detype, Iterative)
        @warn("Beddoes Leishman Iterative implementation isn't prepared yet. - initialize()")
    else #Model is indicial
        if dsmodel.version==2 # Original AeroDyn Implementation
            return updateenvironment_ADO(p, U, aoa)
        elseif dsmodel.version==3 # AeroDyn with Gonzalez modifications. 
            return updateenvironment_ADG(p, U, aoa)
        end
    end
end


#### Inplace Functions
function (model::BeddoesLeishman)(x, xnew, p, t, dt) 
    if isa(model.detype, Functional)
        @warn("Functional implementation not yet prepared.")
    elseif isa(model.detype, Indicial)
        if model.version==1
            @warn("Original indicial Beddoe-Leishman not prepared for use yet.")
        elseif model.version==2
            error("The original AeroDyn implementation doesn't function in-place yet. ")
            ns = numberofstates(model)
            newstates = Array{eltype(p), 1}(undef, ns)
            for i = 1:model.n
                ps = view(p, 18*(i-1)+1:18*i)
                #[c, a, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, xcp]
                c, a, dcndalpha, alpha0, _, _, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, _, U, aoa = ps #Inputs  
                
                xs = view(x, 22*(i-1)+1:22*i)
                xsnew = view(xnew, 22*(i-1)+1:22*i)

                idx = 22*(i-1)+1:22*i
                newstates[idx] = update_states_ADO!(model, xs, c, a, U, dt, aoa, dcndalpha, alpha0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, i)
            end
            return newstates

        elseif model.version==3
            
            for i = 1:model.n
                ps = view(p, 22*(i-1)+1:22*i)
    
                c, a, dcndalpha, alpha0, _, _, A1, A2, A5, b1, b2, b5, Tf0, Tv0, Tp, Tvl, Tsh, Cn1, _, zeta, U, aoa = ps #Inputs  
                
                xs = view(x, 32*(i-1)+1:32*i) #states for the given section 
                xsnew = view(xnew, 32*(i-1)+1:32*i)

                # idx = 32*(i-1)+1:32*i
                update_states_ADG!(model, xsnew, xs, c, a, U, dt, aoa, dcndalpha, alpha0, A1, A2, A5, b1, b2, b5, Tf0, Tv0, Tp, Tvl, Tsh, Cn1, zeta, i)
            end

        elseif model.version==4
            @warn("AeroDyn Beddoe-Leishman with Minema's modifications not prepared for use yet.")
        end
    end
end








