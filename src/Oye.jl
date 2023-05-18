#=
The Oye model, as given by ???. 

Indicial States
1 - f

Environmental Parameters (y)
U, alpha

=#
export Oye

"""
    Oye(detype::DEType, n::Int, airfoils::Array{Airfoil, 1})

The Øye model struct. It stores airfoil data for every section to be simulated. It can be used as a method to return updated states or state rates depending on it's DEType. 

### Inputs
- detype - The type of model it is, Functional(), Iterative(), or Indicial().
- cflag::Int - A flag to apply the separation delay to the coefficient of 1) lift, 2) normal force. 
- version::Int - A flag to say whether to use 1) Hansen 2008, 2) Larsen's Hermite interpolation from Faber's 2018's implementation of the model, or 3) Øye's parabola fit.
- A::Float - Dynamic stall coefficient. 
"""
struct Oye{TI, TF} <: DSModel
    detype::DEType 
    cflag::TI 
    version::TI
    A::TF

    function Oye{TI, TF}(detype::DEType, cflag::TI, version::TI, A::TF) where {TI, TF}
        if cflag>2
            error("Øye: the cflag only accepts flags of 1 or 2.")
        elseif version>3
            error("Øye: the version only accepts values of 1, 2, or 3.")
        end
        return new{TI, TF}(detype, cflag, version, A)
    end

    Oye(detype::DEType, cflag::TI, version::TI, A::TF) where {TI, TF} = Oye{TI, TF}(detype, cflag, version, A)
end

# function Oye(detype::DEType; cflag=1, version=1, A=4.0)
#     return Oye(detype, cflag, version, A)
# end


function numberofstates(dsmodel::Oye)
    return 1
end

function numberofparams(dsmodel::Oye)
    return 2
end

"""
    get_cn(airfoil)

Determines whether the aerodynamic coefficient should be with respect to the lift force or normal force based off the user's choice.

**Arguments**
- airfoil::Airfoil: The airfoil being evaluated.
"""
function get_cn(airfoil, alpha)
    if airfoil.model.cflag == 2
        return airfoil.cn(alpha) #Static normal force
    else
        return airfoil.cl(alpha)
    end
end

"""
    get_dcndalpha(airfoil)

Determines whether the derivative of the aerodynamic coefficient should be with respect to the lift force or normal force based off the user's choice.

**Arguments**
- airfoil::Airfoil: The airfoil being evaluated.
"""
function get_dcndalpha(airfoil)
    if airfoil.model.cflag == 2
        return airfoil.dcndalpha #Static normal force
    else
        return airfoil.dcldalpha
    end
end

# function (dsmodel::Oye)(x, p, t, dt)
#     if isa(dsmodel.detype, Functional)
#         error("The state space Oye model is not setup yet.")
#     elseif isa(dsmodel.detype, Iterative)
#         error("The iterative Oye model is not set up yet.")
#     else #The model is indicial
#         nst = numberofstates_total(dsmodel)
#         ns = numberofstates(dsmodel)
#         np = numberofparams(dsmodel)
#         newstates = Array{eltype(p), 1}(undef, nst)
#         for i = 1:dsmodel.n
#             ps = view(p, np*(i-1)+1:np*i)
        
#             c, dcndalpha, alpha0, alphasep, A, U, aoa = ps #Inputs 
            
#             xs = view(x, ns*(i-1)+1:ns*i) #Nodal states. 
            
#             idx = ns*(i-1)+1:ns*i
#             newstates[idx] = update_states(dsmodel, xs, U, aoa, dt, c, dcndalpha, alpha0, alphasep, A, i)
#         end
#         return newstates
#     end
# end

"""
    initialize(dsmodel::Oye, airfoil::Airfoil, tvec, y)

Finds the initial state value for an airfoil that will be used in the indicial solve method.

**Arguments**
- dsmodel::DSModel: The model being used.
- airfoil::Airfoil: The airfoil being evaluated.
- tvec::StepRangeLen: The range of time values that are used for calculations.
- y::Vector{TF}: A vector containing the parameters for an airfoil (such as the angle of attack at a given time).
"""
function initialize(dsmodel::Oye, airfoil::Airfoil, tvec, y)
    if isa(dsmodel.detype, Functional)
        @warn("Oye Functional implementation isn't prepared yet. - initialize()")
    elseif isa(dsmodel.detype, Iterative)
        @warn("Oye Iterative implementation isn't prepared yet. - initialize()")
    else #Model is indicial

        _, _, alpha, _ = y

        fst = separationpoint(airfoil, alpha)

        if fst>1
            fst = 1.0
        end
        
        states = [fst]

        loads = zeros(3)
        get_loads!(dsmodel::Oye, airfoil, states, loads, y)
        
        return states, loads, y
    end
end






######## Indicial code
function update_states(dsmodel::Oye, airfoil::Airfoil, oldstate, y, dt)
    newstate = zero(oldstate)
    update_states!(dsmodel, airfoil, oldstate, newstate, y, dt)
    return newstate
end

#=
In-place version. 
=#
"""
    update_states!(dsmodel::Oye, airfoil::Airfoil, oldstate, newstate, y, dt)

Finds the new state values for an arifoil at a given time value using the indicial solve method.

**Arguments**
- dsmodel::DSModel: The model being used.
- airfoil::Airfoil: The airfoil being evaluated.
- oldstate::Vector{TF}: The old/previous states for the evaluated airfoil.
- newstate::Vector{TF}: The new state values that are solved for from this function.
- y::Vector{TF}: A vector containing the parameters for an airfoil (such as the inflow velocity).
- dt::TF: The value of the magnitude of the time step for an airfoil.
"""
function update_states!(dsmodel::Oye, airfoil::Airfoil, oldstate, newstate, y, dt)
    ### Unpack 
    U, _, alpha, _ = y
    fold = oldstate[1]


    fst = separationpoint(airfoil, alpha) #Current static degree of attachment (separation point)  

    tau = dsmodel.A*airfoil.c/U #Time constant - From Hansen's 2008 paper (all other constants will need to be converted to Hansen's format)

    f = fst + (fold - fst)*exp(-dt/tau) #Delay on separation point 

    if f>1
        f = 1.0
    end

    newstate[1] = f
end

"""
    get_loads(dsmodel::Oye, airfoil::Airfoil, states, y)

Obtains the coefficients of lift, drag, and moments for a given airfoil using the indicial solve method.

**Arguments**
- dsmodel::DSModel: The model being used.
- airfoil::Airfoil: The airfoil being evaluated.
- states::Vector{TF}: The states at a given time value.
- y::Vector{TF}: A vector containg the parameters for an airfoil (such as its angle of attack).
"""
function get_loads(dsmodel::Oye, airfoil::Airfoil, states, y)
    ### Unpack
    f = states[1]

    _, _, alpha, _ = y # U, alpha = y


    dcndalpha = get_dcndalpha(airfoil)

    alpha0 = airfoil.alpha0

    cn_inv = dcndalpha*(alpha-alpha0) #Hansen 2004 EQ 19.

    if dsmodel.version==2 #Larsen (Faber) #TODO: I wonder if there is a good way to use multiple dispatch on this. 
        cn_fs = cl_fullysep_faber(airfoil, alpha) #this was originally cl_fullysep_larsen, but it needs to be cl_fullysep_faber
    elseif dsmodel.version==1 #Hansen
        cn_fs = cl_fullysep_hansen(airfoil, alpha)
    else #Oye Parabola Fit
        cn_fs = cl_fullysep_oye(airfoil, alpha)
    end


    #Todo: Calculate the Cc, Cl, Cd, and potentially the Cm if possible. 
    if dsmodel.cflag==2 #delay applied to normal and tangential loads
        #Rotate loads
    else 
        Cl = f*cn_inv + (1-f)*cn_fs #Hansen 2004 EQ 17. Coefficient of lift
        Cd = zero(Cl) #Coefficient of drag
        Cm = zero(Cl) #Coefficient of moment
    end

    return Cl, Cd, Cm
end

function get_loads!(dsmodel::Oye, airfoil::Airfoil, states, loads, y)
    loads .= get_loads(dsmodel, airfoil, states, y)
end

"""
    cl_fullysep_hansen(airfoil, alpha)

Finds the fully separated coefficient of lift values using the method found in Hansen's papers.

**Arguments**
- airfoil::Airfoil: The airfoil being evaluated
- alpha::TF: The angle of attack value being passed in.
"""
function cl_fullysep_hansen(airfoil, alpha) #Todo: Move to Hansen's file
    cn = get_cn(airfoil, alpha) #Static lift
    dcndalpha = get_dcndalpha(airfoil)

    cn_inv = dcndalpha*(alpha-airfoil.alpha0) #Hansen 2004 EQ 19, inviscid lift
    
    # fst = (2*sqrt(abs(cn/cn_inv)) - 1)^2 #Hansen 2004 EQ 15.
    fst = separationpoint(airfoil, alpha)
    cn_fs = (cn-(cn_inv*fst))/(1-fst) #Hansen 2004 EQ 18, fully separated lift

    return cn_fs
end

#=
Larsen's fully separated lift coefficient from Faber's 2018 masters thesis. 
=#
"""
    cl_fullysep_faber(airfoil, alpha)

Finds the fully separated coefficient of lift values using the Hermite interpolation method found in Faber's paper.

**Arguments**
- airfoil::Airfoil: The airfoil being evaluated.
- alpha::TF: The angle of attack value being passed in.
"""
function cl_fullysep_faber(airfoil, alpha) #Todo: Move to Larsen's file
    alpha0 = airfoil.alpha0
    alpha_sep = airfoil.alphasep[2]

    cn = get_cn(airfoil, alpha)
    cn_sep = get_cn(airfoil, alpha_sep)
    dcndalpha = get_dcndalpha(airfoil)

    if alpha0 < alpha < alpha_sep
        #Hermite Interpolation as in Faber 2018
        #TODO: Needs to be extended to work for angles less than alpha0, it should have a reflection to what is done here. 
        t0 = (alpha - alpha0)/(alpha_sep - alpha0) #Faber 2018 EQ A.1a
        #As alpha -> alpha_sep this will diverge. So... what keeps it from diverging? 
        t1 = (alpha - alpha_sep)/(alpha_sep - alpha0) #Faber 2018 EQ A.1b
        term1 = (alpha_sep - alpha0)*dcndalpha*(1 + t0*(7*t1/6 - 1))/2
        term2 = cn_sep*t0*(1-( 2*t1))
        clfs = t0*(term1+term2) #Faber 2018 EQ A.4

        return clfs
    else
        return cn
    end
end

"""
    cl_fullysep_oye(airfoil, alpha)

Finds the fully separated coefficient of lift values using Øye's orginal method.

**Arguments**
- airfoil::Airfoil: The airfoil being evaluated.
- alpha::TF: The angle of attack value being passed in.
"""
function cl_fullysep_oye(airfoil, alpha)
    
    a0 = airfoil.alpha0
    a_sep = airfoil.alphasep[2]
    Cl_sep = airfoil.cl(a_sep)
    dclda = airfoil.dcldalpha * 0.5

    
    p1 = (a0*dclda + Cl_sep - dclda*a_sep)/((a_sep - a0)^2)
    p2 = (-dclda*a0^2 + dclda*a_sep^2 - 2*Cl_sep*a0)/((a_sep - a0)^2)
    p3 = (dclda*a0^2*a_sep - Cl_sep*a0^2 - dclda*a_sep^2*a0 + 2*a0^2*Cl_sep)/((a_sep - a0)^2)


    if a0 <= alpha <= a_sep
        Cl_fullysep = p1*alpha^2 + p2*alpha + p3
    else 
        Cl_fullysep = airfoil.cl(alpha)
    end


    return Cl_fullysep
end




"""
    state_rates!(model::Oye, airfoil::Airfoil, dx, x, y, t)

Calculate the state rates of the Øye model. 

**Arguments**
- model::DSModel: The model to be used.
- airfoil::Airfoil: The airfoil to be used.
- dx::Vector: The vector to be filled with the state rates.
- x::Vector: The current state.
- y::Vector: The current input.
- t::Float64: The current time.
"""
function state_rates!(model::Oye, airfoil::Airfoil, dx, x, y, t)
    
    ### Evaluate environmental functions
    U, _, alpha, _ = evaluate_environment(y, t)
    

    ### Prepare time constant
    A = model.A
    c = airfoil.c
    T_f = U/(A*c)

    ### Fetch static separation point
    fst = separationpoint(airfoil, alpha)
    # @show fst

    ### Calculate state rate
    dx[1] = -T_f*(x[1]-fst) 
end

export parsesolution

"""
    parsesolution(model::Oye, airfoils::AbstractVector{<:Airfoil}, sol, y)

Calculates the dynamic lift coefficient values for simulated airfoils and gives their corresponding angle of attacks.

**Arguments**
- model::DSModel: The model to be used.
- airfoils::AbstractVector{<:Airfoil}: The vector of airfoils being evaluated.
- sol::ODESolution: The solutions to the state rate problem found using DifferentialEquations.jl.
- y::Vector{TF}: The vector containing the parameter values for the ODE problem.

**Outputs**
- Lift_aoa_Matrix::Matrix{TF}: The matrix containing the dynamic lift coefficient values and the corresponding angle of attacks for each airfoil.
"""
function parsesolution(model::Oye, airfoils::AbstractVector{<:Airfoil}, sol::ODESolution, y)
    f = Array(sol) 
    tvec = sol.t 


    #preallocates a matrix to be filled with the dynamic lift and angle of attack values for each airfoil evaluated
    Lift_aoa_Matrix = zeros(2*length(airfoils) , length(tvec))  #needs to be changed for a list of airfoils


    for w in 1:length(airfoils) #this needs to be changed to loop through all of the airfoils
        airfoil = airfoils[w]

        alpha0 = airfoil.alpha0 

        alphavec = y[3+(w-1)*4]

        if airfoil.model.cflag == 1 #checks to see if the user chose the coefficient of lift flag
            cl_sep = airfoil.cl(airfoil.alphasep[2])
        else
            cl_sep = airfoil.cn(airfoil.alphasep[2]) #find the coefficient of normal force at the fully seperated angle of attack 
        end



        for i in 1:length(tvec) 


            Lift_aoa_Matrix[2*w-1, i] = alphavec(tvec[i]) 
            C_inv = airfoil.dcldalpha*(alphavec(tvec[i]) - alpha0) 



            if airfoil.model.cflag == 1 #checks to see if the user used the coefficient of lift flag
                cl = airfoil.cl(alphavec(tvec[i])) 
            else
                cl = airfoil.cn(alphavec(tvec[i])) #finds what the coefficient of static normal force value is at a given angle of attack
            end




            if airfoil.model.version == 1 #checks to see if the Hansen method is chosen for this Oye solve
                fst = (2*sqrt(abs(cl/C_inv))-1)^2 
                C_fs = (cl - (C_inv*fst))/(1-fst) #finds the fully separated coefficient of lift value 
            elseif airfoil.model.version == 2 #for the faber method of fully separated lift
                C_fs = cl_fullysep_faber(airfoil, alphavec(tvec[i])) #finds the fully separated coefficient of lift value for the Faber method
            else #for the Øye way of finding fully separated lift
                C_fs = cl_fullysep_oye(airfoil, alphavec(tvec[i]))
            end



            

            C_L_dyn = f[w,i]*C_inv+(1-f[w,i])*C_fs 
            Lift_aoa_Matrix[2*w,i] = C_L_dyn 
        end

    end

    return Lift_aoa_Matrix
end


