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
struct BeddoesLeishman{TI, TF} <: DSModel
    detype::DEType 
    version::TI #Which version of the indicial implementation.  
    A #::TFV #The dynamic pressure response coefficients
    b #::TFV #The secondary dynamic pressure coefficients
    T #::TFV #The time constants
    Cn1::TF #Separation normal force coefficient
    Cd0::TF #Zero lift drag (Viscous drag)
    Cm0::TF #Zero lift moment
    eta::TF #The recovery efficiency
    zeta::TF #I forget which efficiency this is. 
    a::TF # Speed of sound

    function BeddoesLeishman{TI, TF}(detype::DEType, version::TI, A, b, T, Cn1::TF, Cd0::TF, Cm0::TF, eta::TF, zeta::TF, a::TF) where {TI, TF}
        if version>4
            error("BeddoesLeishamn: the version only accepts whole integers between 1-4.")
        end
        return new{TI, TF}(detype, version, A, b, T, Cn1, Cd0, Cm0, eta, zeta, a)
    end

    BeddoesLeishman(detype::DEType, version::TI, A, b, T, Cn1::TF, Cd0::TF, Cm0::TF, eta::TF, zeta::TF, a::TF) where {TI, TF} = BeddoesLeishman{TI, TF}(detype, version, A, b, T, Cn1, Cd0, Cm0, eta, zeta, a)
end

function BeddoesLeishman(detype::DEType, version::Int, polar, alpha0, alpha1, A, b, T; eta=1.0, zeta=0.5, a=335.0, interp=Akima)
    clfit = interp(polar[:,1], polar[:,2])
    cdfit = interp(polar[:,1], polar[:,3])
    cmfit = interp(polar[:,1], polar[:,4])

    Cd0 = cdfit(alpha0)
    Cm0 = cmfit(alpha0)
    Cn1 = clfit(alpha1)*cos(alpha1) + (cdfit(alpha1) - Cd0)*sin(alpha1)

    return BeddoesLeishman(detype, version, A, b, T, Cn1, Cd0, Cm0, eta, zeta, a)
end



function getloads(dsmodel::BeddoesLeishman, states, p, airfoil)
    loads = zeros(3)
    getloads!(dsmodel, airfoil, states, loads, y)
end

function getloads!(dsmodel::BeddoesLeishman, airfoil::Airfoil, states, loads, y)
    if isa(dsmodel.detype, Functional)
        @warn("Functional implementation not yet prepared.")
    elseif isa(dsmodel.detype, Indicial)
        if dsmodel.version==1
            @warn("Original indicial Beddoe-Leishman not prepared for use yet.")
        elseif dsmodel.version==2
            return getloads_BLA(dsmodel, states, p, airfoil)
        elseif dsmodel.version==3
            # @warn("AeroDyn Beddoe-Leishman with Gonzalez's modifications not prepared for use yet.")
            # return getloads_BLAG(dsmodel, states, p, airfoil)
            BLADG_coefficients!(airfoil::Airfoil, loads, states, y)
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

# function initialize(dsmodel::Oye, airfoil::Airfoil, tvec, y)
function initialize(dsmodel::BeddoesLeishman, airfoil::Airfoil, tvec, y)
    if isa(dsmodel.detype, Functional)
        @warn("Beddoes Leishman Functional implementation isn't prepared yet. - initialize()")
    elseif isa(dsmodel.detype, Iterative)
        @warn("Beddoes Leishman Iterative implementation isn't prepared yet. - initialize()")
    else #Model is indicial
        if dsmodel.version==2 # Original AeroDyn Implementation
            return initialize_ADO(Uvec, aoavec, tvec, airfoil, c, a)
        elseif dsmodel.version==3 # AeroDyn with Gonzalez modifications. 
            return initialize_ADG(airfoil, tvec, y)
        end
    end
end



function update_states(dsmodel::BeddoesLeishman, airfoil::Airfoil, x, y, dt) 
    ns = numberofstates(model)
    newstates = Array{eltype(p), 1}(undef, ns)

    update_states!(dsmodel, airfoil, x, newstates, y, dt)
            
    return newstates
end

#### Inplace Functions
function update_states!(model::BeddoesLeishman, airfoil::Airfoil, x, xnew, y, dt) 
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

                update_states_ADG!(airfoil, x, xnew, y, dt)

        elseif model.version==4
            @warn("AeroDyn Beddoe-Leishman with Minema's modifications not prepared for use yet.")
        end
    end
end








