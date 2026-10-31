

export BeddoesLeishman

"""
    BeddoesLeishman(detype, version, polar, alpha0, alpha1, A, b, T;
                    eta=1.0, zeta=0.5, a=335.0, interp=Akima)

Model constants for one airfoil. The Gonzalez example uses OpenFASTTools to
construct a Discrete() version 3 model and its Airfoil together. It does not
use the obsolete vector-of-airfoils model constructor.

### Inputs
- detype - The type of model it is, Continuous(), Discrete(), or Discrete().
- polar - Columns alpha (radians), Cl, Cd, Cm.
- alpha0, alpha1 - Zero-lift and positive separation angles (radians).
- A, b, T - Response amplitudes, decay constants and model time constants.
- version - Which version of the Discrete implementation. 1) original. 2) AeroDyn original. 3) AeroDyn Gonzalez. 4) AeroDyn Minema
"""
struct BeddoesLeishman{TI, TF, TFV} <: DSModel
    detype::DEType 
    version::TI #Which version of the Discrete implementation.  
    A::TFV #The dynamic pressure response coefficients #Todo: I didn't write down why I changed these from TFV. 
    b::TFV #The secondary dynamic pressure coefficients
    T::TFV #The time constants
    Cn1::TF #Separation normal force coefficient
    Cd0::TF #Zero lift drag (Viscous drag)
    Cm0::TF #Zero lift moment
    eta::TF #The recovery efficiency
    zeta::TF #I forget which efficiency this is. 
    a::TF # Speed of sound

    function BeddoesLeishman{TI, TF, TFV}(detype::DEType, version::TI, A::TFV, b::TFV, T::TFV, Cn1::TF, Cd0::TF, Cm0::TF, eta::TF, zeta::TF, a::TF) where {TI, TF, TFV}
        if !(version in 1:4)
            error("BeddoesLeishamn: the version only accepts whole integers between 1-4.")
        end
        isfinite(a) && a>0 || throw(ArgumentError("Sound speed must be finite and positive."))
        return new{TI, TF, TFV}(detype, version, A, b, T, Cn1, Cd0, Cm0, eta, zeta, a)
    end

    BeddoesLeishman(detype::DEType, version::TI, A::TFV, b::TFV, T::TFV, Cn1::TF, Cd0::TF, Cm0::TF, eta::TF, zeta::TF, a::TF) where {TI, TF, TFV} = BeddoesLeishman{TI, TF, TFV}(detype, version, A, b, T, Cn1, Cd0, Cm0, eta, zeta, a)
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



function get_loads!(model::BeddoesLeishman, airfoil::Airfoil, states, loads, y, p)
    require_adg(model)
    BLADG_coefficients!(airfoil,loads,states,y,p)
    return loads
end

function get_loads(model::BeddoesLeishman, airfoil::Airfoil, states, y, p)
    require_adg(model)
    return BLADG_coefficients(airfoil,states,y,p)
end

function numberofstates(dsmodel::BeddoesLeishman) #TODO: This probably need to be augmented to check if the model is a Continuous, an Discrete, or an Discrete. 
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
    else
        dsv = dsmodel.version
        error("$dsv isn't a valid version model number for the Beddoes-Leishman models.")
        return 0
    end
end

function numberofloads(dsmodel::BeddoesLeishman) #TODO: This probably need to be augmented to check if the model is a Continuous, an Discrete, or an Discrete.
    if dsmodel.version==1
        @warn("The orginal Beddoes-Leishman model is not yet prepared.")
        return 0
    elseif dsmodel.version==2
        return 5
    elseif dsmodel.version==3
        # @warn("The Gonzalez Beddoes-Leishman model is not yet prepared.")
        return 3
    elseif dsmodel.version==4
        @warn("The Minema Beddoes-Leishman model is not yet prepared.")
        return 5
    else
        dsv = dsmodel.version
        error("$dsv isn't a valid version model number for the Beddoes-Leishman models.")
        return 0
    end
end

function numberofparams(dsmodel::BeddoesLeishman) #TODO: This probably need to be augmented to check if the model is a Continuous, an Discrete, or an Discrete.

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
function initialize(model::BeddoesLeishman, airfoil::Airfoil, tvec, y, p)
    require_adg(model)
    return initialize_ADG(airfoil,tvec,y,p)
end

function update_states(dsmodel::BeddoesLeishman, airfoil::Airfoil, x, y, p, dt) 
    require_adg(dsmodel)
    ns = numberofstates(dsmodel)
    newstates = similar(x,promote_type(eltype(x),eltype(airfoil.polar)),ns)

    update_states!(dsmodel, airfoil, x, newstates, y, p, dt)
            
    return newstates
end

#### Inplace Functions
function update_states!(model::BeddoesLeishman, airfoil::Airfoil, x, xnew, y, p, dt)
    require_adg(model)
    update_states_ADG!(airfoil,x,xnew,y,p,dt)
    return xnew
end

function require_adg(model::BeddoesLeishman)
    model.detype isa Discrete && model.version==3 || throw(ArgumentError("Only Discrete() BeddoesLeishman version 3 (Gonzalez) is implemented. Other variants are unfinished."))
    return nothing
end

function state_rates!(model::BeddoesLeishman, airfoil::Airfoil, dx, x, y, p, t)
    throw(ArgumentError("Continuous BeddoesLeishman is unfinished; use Discrete() version 3."))
end
