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
- version::Int - A flag to say whether to use 1) Hansen 2008, or 2) Larsen's Hermite interpolation from Faber's 2018's implementation of the model.
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
        elseif version>2
            error("Øye: the version only accepts values of 1 or 2.")
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

function get_cn(airfoil, alpha)
    if airfoil.model.cflag == 2
        return airfoil.cn(alpha) #Static normal force
    else
        return airfoil.cl(alpha)
    end
end

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

function get_loads(dsmodel::Oye, airfoil::Airfoil, states, y)
    ### Unpack
    f = states[1]

    _, _, alpha, _ = y # U, alpha = y


    dcndalpha = get_dcndalpha(airfoil)

    alpha0 = airfoil.alpha0

    cn_inv = dcndalpha*(alpha-alpha0) #Hansen 2004 EQ 19.

    if dsmodel.version==2 #Larsen (Faber) #TODO: I wonder if there is a good way to use multiple dispatch on this. 
        cn_fs = cl_fullysep_faber(airfoil, alpha) #this was originally cl_fullysep_larsen, but it needs to be cl_fullysep_faber
    else #Hansen
        cn_fs = cl_fullysep_hansen(airfoil, alpha)
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


function parsesolution(model::Oye, airfoil::Airfoil, sol, y)
    _, _, alphavec, _ = y
    f = Array(sol) 
    tvec = sol.t 


    #preallocates a matrix to be filled with the dynamic lift and angle of attack values for each airfoil evaluated
    Lift_aoa_Matrix = zeros(2*1 , length(tvec))  #needs to be changed for a list of airfoils



    for w in 1:1 #this needs to be changed to loop through all of the airfoils

        alpha0 = airfoil.alpha0 



        if airfoil.model.cflag == 1 #checks to see if the user chose the coefficient of lift flag
            cl_sep = airfoil.cl(airfoil.alphasep[2])
        else
            cl_sep = airfoil.cn(airfoil.alphasep[2]) #find the coefficient of normal force at the fully seperated angle of attack 
        end



        for i in 1:eachindex(tvec) 


            Lift_aoa_Matrix[w, i] = alphavec(tvec[i]) 
            C_inv = airfoil.dcldalpha*(alphavec(tvec[i]) - alpha0) 



            if airfoil.model.cflag == 1 #checks to see if the user used the coefficient of lift flag
                cl = airfoil.cl(alphavec(tvec[i])) 
            else
                cl = airfoil.cn(alphavec(tvec[i])) #finds what the coefficient of static normal force value is at a given angle of attack
            end




            if airfoil.model.version == 1 #checks to see if the Hansen method is chosen for this Oye solve
                fst = (2*sqrt(abs(cl/C_inv))-1)^2 
                C_fs = (cl - (C_inv*fst))/(1-fst) #finds the fully separated coefficient of lift value 
            elseif airfoil.model.version == 2
                C_fs = cl_fullysep_faber(airfoil, alphavec(tvec[i])) #finds the fully separated coefficient of lift value for the Faber method
            end




            C_L_dyn = f[w,i]*C_inv+(1-f[w,i])*C_fs 
            Lift_aoa_Matrix[w+1,i] = C_L_dyn 
        end

    end

    return Lift_aoa_Matrix 

end


