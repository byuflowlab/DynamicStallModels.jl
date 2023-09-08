export Riso

"""
    Riso(detype::DEType, A::Vector, b::Vector, T::vector)

The Riso model struct. This struct holds information pertaining to the solution approach and the coefficients that are used in this model.

### Inputs
- detype - The type of solution approach used, Indicial(), or Functional()
- A::Vector - A vector containing the coefficients for A1 and A2
- b::Vector - A vector containing the coefficients for b1 and b2
- T::Vector - A vector containing the coefficients for Tp and Tf
"""
struct Riso<: DSModel
    detype::DEType
    A
    b
    T
end


function numberofstates(dsmodel::Riso)
    if isa(dsmodel.detype, Indicial) #Indicial model contains extra states to account for previous time step parameter values
        return 8
    else
        return 4 #This is for the Functional method
    end
end

function numberofparams(dsmodel::Riso)
    return 6
end



"""
    initialize(dsmodel::Riso, airfoil::Airfoil, tvec, y)

Finds the initial state value for an airfoil that will be used in the indicial solve method.

**Arguments**
- dsmodel::DSModel: The model being used.
- airfoil::Airfoil: The airfoil being evaluated.
- tvec::StepRangeLen: The range of time values that are used for calculations.
- y::Vector{TF}: A vector containing the parameters for an airfoil (such as the angle of attack at a given time).
"""
function initialize(dsmodel::Riso, airfoil::Airfoil, tvec, y)
    if isa(dsmodel.detype, Indicial)
        U, Udot, alpha, alphadot = y
        states = [0.0, 0.0, 0.0, 0.0, U, Udot, alpha, alphadot] #The initial parameters need to become states for the Indicial formulation

        loads = zeros(3)
        get_loads!(dsmodel::Riso, airfoil::Airfoil, states, loads, y)

        return states, loads, y
    end
end


function update_states(model::Riso, airfoil::Airfoil, oldstate, y, dt)
    newstate = zero(oldstate)
    update_states!(model::Riso, airfoil::Airfoil, oldstate, newstate, y, dt)
    return newstate
end


"""
    update_states!(model::Riso, airfoil::Airfoil, oldstate, newstate, y, dt)

    Solves the new state values of after a given time step during the indicial solve.

**Arguments**
- model::Riso - The Dynamic Stall model being used.
- airfoil::Airfoil - The airfoil being evaluated.
- oldstate - A vector containing the previous time step states and also values like angle of attack, inflow velocity, etc.
- newstate - A vector that will be filled with the newly calculated state values at the end of the function.
- y - A vector containing parameter values like the angle of attack and inflow velocity at the new time step.
- dt - The magnitude of the length of time in between time steps.
"""
function update_states!(model::Riso, airfoil::Airfoil, oldstate, newstate, y, dt)
    U_new, Udot_new, alpha_new, alphadot_new = y #Riso's method requires parameter values from the current and previous time steps
    x1_old , x2_old, x3_old, x4_old, U_old, Udot_old, alpha_old, alphadot_old = oldstate 


    A1 = model.A[1] #Riso's model carries these coefficients when the model is created
    A2 = model.A[2]
    b1 = model.b[1]
    b2 = model.b[2]


    c = airfoil.c
    dcldalpha = airfoil.dcldalpha
    alpha0 = airfoil.alpha0
    Tu = c/(2*U_new)
    Tp = Tu*model.T[1] #T[1] and T[2] are just given coefficient values that are said to be good, but these can change as well
    Tf = Tu*model.T[2]


    ae_old = alpha_old*(1-A1-A2) + x1_old + x2_old


    P1 = b1*(U_new + U_old)/c + (Udot_new + Udot_old)/(U_new + U_old)
    P2 = b2*(U_new + U_old)/c + (Udot_new + Udot_old)/(U_new + U_old)
    P3 = 1/Tp
    P4 = 1/Tf


    Q1 = (U_old*alpha_old + U_new*alpha_new)*(b1*A1)/(c)
    Q2 = (U_old*alpha_old + U_new*alpha_new)*(b2*A2)/(c)


    x1_new = exp(-1*P1*dt)*x1_old + (Q1/P1)*(1 - exp(-1*P1*dt))
    x2_new = exp(-1*P2*dt)*x2_old + (Q2/P2)*(1 - exp(-1*P2*dt))


    ae_new = alpha_new*(1-A1-A2) + x1_new + x2_new


    Q3 = (0.5/Tp)*(dcldalpha*(ae_old + ae_new - 2*alpha0) + (pi*c*0.5)*((alphadot_old/U_old + alphadot_new/U_new))) 
    x3_new = exp(-1*P3*dt)*x3_old + (Q3/P3)*(1 - exp(-1*P3*dt))


    Q4 = (0.5/Tf)*(separationpoint(airfoil, x3_old/(dcldalpha+alpha0)) + separationpoint(airfoil, x3_new/(dcldalpha+alpha0)))


    x4_new = exp(-1*P4*dt)*x4_old + (Q4/P4)*(1 - exp(-1*P4*dt))


    newstate[1] = x1_new
    newstate[2] = x2_new
    newstate[3] = x3_new
    newstate[4] = x4_new
    newstate[5] = U_new #The current parameters are pushed into the newstate vector, so they can become old state values (previous time step)
    newstate[6] = Udot_new
    newstate[7] = alpha_new
    newstate[8] = alphadot_new
end

"""
    get_loads(model::Riso, airfoil::Airfoil, states, y)

Calculates the coefficients of lift, drag, and moment at a given time step for the Indicial method.

**Arguments**
- model::Riso - The Dynamic Stall model used in the evaluation.
- airfoil::Airfoil - The specific airfoil being evaluated.
- states - The states at a given time step.
- y - The parameters of a given time step: such as inflow velocity and angle of attack.
"""
function get_loads(model::Riso, airfoil::Airfoil, states, y)
    x1, x2, _, x4, _, _, _, _ = states #unpacking

    U, _, alpha, alphadot = y #only current parameter values are needed to find the loads

    alpha0 = airfoil.alpha0
    c = airfoil.c
    dcldalpha = airfoil.dcldalpha
    CD0 = airfoil.cd(alpha0) #CD0 and CM0 are the drag and moment coefficient values at the zero lift angle of attack
    CM0 = airfoil.cm(alpha0)

    A1 = model.A[1]
    A2 = model.A[2]


    ae = alpha*(1 - A1 - A2) + x1 + x2
    Tu = c/(U*2)
    CL_fs = hansen_fully_sep(airfoil , ae) #Fully separated lift value at the effective angle of attack
    CD_st = airfoil.cd(ae) #The static drag and moment coefficients at the effective angle of attack
    CM_st = airfoil.cm(ae)
    ast_1 = (airfoil.cm(x4 - CM0)/airfoil.cl(x4)) #ast_1 and ast_2 are passed in to find the dynamic coefficients
    ast_2 = (airfoil.cm(separationpoint(airfoil, ae)) - CM0)/airfoil.cl(separationpoint(airfoil, ae)) 

    CL_Dyn = dcldalpha*(ae - alpha0)*x4 + CL_fs*(1 - x4) + pi*Tu*alphadot
    CD_Dyn = CD_st + (alpha - ae)*CL_Dyn + (CD_st - CD0)*((sqrt(separationpoint(airfoil, ae)) - sqrt(x4))/2 - (separationpoint(airfoil,ae) - x4)/4)
    CM_Dyn = CM_st + CL_Dyn*(ast_1 - ast_2) - (pi/2)*Tu*alphadot

    return CL_Dyn, CD_Dyn, CM_Dyn

end

function get_loads!(model::Riso, airfoil::Airfoil, states, loads, y)
    loads .= get_loads(model::Riso, airfoil::Airfoil, states, y)
end


"""
    state_rates!(model::Riso, airfoil::Airfoil, dx, x, y, t)

Solves for the state values using the Riso method for a specific airfoil over a range of time values.

**Arguments**
- model::Riso - The Dynamic Stall model being used to evaluate.
- airfoil::Airfoil - The specific airfoil being evaluated.
- dx - A vector containing the state rate values.
- x - A vector containing the state values.
- y - A vector that contains functions used to calculate the parameters at a given time.
- t - The time value that the solve is on.
"""
function state_rates!(model::Riso , airfoil::Airfoil , dx, x, y, t)
    
    U, Udot, alpha, alphadot = evaluate_environment(y,t) #Finds the parameter values at a given time


    A1 = model.A[1]
    A2 = model.A[2]
    b1 = model.b[1]
    b2 = model.b[2]


    c = airfoil.c
    dcldalpha = airfoil.dcldalpha
    alpha0 = airfoil.alpha0
    Tu = c/(2*U)
    Tp = Tu*model.T[1]
    Tf = Tu*model.T[2]

    ae = alpha*(1-A1-A2) + x[1] + x[2] #Effective angle of attack

    fst_val = separationpoint(airfoil, (alpha0 + x[3]/dcldalpha))

    dx[1] = b1*A1*alpha/Tu - x[1]*(b1+c*Udot/(2*U^2))/Tu #State rate equations for Riso
    dx[2] = b2*A2*alpha/Tu - x[2]*(b2 + c*Udot/(2*U^2))/Tu
    dx[3] = (dcldalpha*(ae-alpha0) + pi*Tu*alphadot)/Tp - x[3]/Tp
    dx[4] = fst_val/Tf - x[4]/Tf
end


"""
    hansen_fully_sep(airfoil, alpha)

Calculates the fully separated lift value for the Riso method.

**Arguments**
- airfoil::Airfoil - The specific airfoil being evaluated.
- alpha - The angle of attack being evaluated at.
"""
function hansen_fully_sep(airfoil , alpha)
    CL_st = airfoil.cl(alpha) #Finds the static lift coefficient values at the desired angle of attack
    dcldalpha = airfoil.dcldalpha
    alpha0 = airfoil.alpha0
    delta_alpha = (alpha-alpha0) #Used to multiply with the linear slope for inviscid lift

    fst = (2*sqrt(CL_st/(dcldalpha*(delta_alpha))) - 1)^2 #This equation is needed for the fully seperated lift equation to not give undefined results

    CL_fs = (CL_st - dcldalpha*delta_alpha*fst)/(1-fst)

    return CL_fs
end

export parsesolution


"""
    parsesolution(model::Riso, airfoils::AbstractVector{<:Airfoil}, sol::ODESolution, y)

Used to find the load values for the Riso model using the Functional method.

**Arguments**
- model::Riso - The Dynamic Stall model being used.
- airfoils::AbstractVector{<:Airfoil} - A vector containing all of the airfoils being evaluated.
- sol::ODESolution - The solution that came from the ODE solve.
- y - A vector containing all of the functions that solve for the parameters of each airfoil.

**Outputs**
- coefficient_matrix::Matrix{TF} - The matrix containing the dynamic lift coefficient values and the corresponding angle of attacks for each airfoil.
"""
function parsesolution(model::Riso, airfoils::AbstractVector{<:Airfoil}, sol::ODESolution, y)
    f = Array(sol) #Changes the solution from the ODE problem into a useable matrix
    tvec = sol.t #Creates a vector of the evaluated times from the ODE problem


    coefficient_matrix = zeros(4*length(airfoils) , length(tvec)) #generates a matrix that is big enough for each airfoil and their corresponding coefficients
    
    for w in 1:length(airfoils) #Loops through each airfoil 

        airfoil = airfoils[w]

        alpha0 = airfoil.alpha0 
        c = airfoil.c
        U = y[1+(w-1)*4] #Obtains the correct parameter values for the corresponding airfoil
        alphadot = y[4+(w-1)*4]
        alpha = y[3+(w-1)*4]
        dcldalpha = airfoil.dcldalpha
        CD0 = airfoil.cd(alpha0)
        CM0 = airfoil.cm(alpha0)

        A1 = model.A[1]
        A2 = model.A[2]


        for i in 1:length(tvec)
           
            ae = alpha(tvec[i])*(1 - A1 - A2) + f[1+(w-1)*4, i] + f[2+(w-1)*4, i]
            Tu = c/(U(tvec[i])*2)
            CL_fs = hansen_fully_sep(airfoil , ae)
            CD_st = airfoil.cd(ae)
            CM_st = airfoil.cm(ae)
            ast_1 = (airfoil.cm(f[4+(w-1)*4, i]) - CM0)/airfoil.cl(f[4+(w-1)*4, i]) 
            ast_2 = (airfoil.cm(separationpoint(airfoil, ae)) - CM0)/airfoil.cl(separationpoint(airfoil, ae)) 



            CL_Dyn = dcldalpha*(ae - alpha0)*f[4+(w-1)*4, i] + CL_fs*(1 - f[4+(w-1)*4, i]) + pi*Tu*alphadot(tvec[i]) #Dynamic coefficients are solved for in this block
            CD_Dyn = CD_st + (alpha(tvec[i]) - ae)*CL_Dyn + (CD_st - CD0)*((sqrt(separationpoint(airfoil, ae)) - sqrt(f[4+(w-1)*4, i]))/2 - (separationpoint(airfoil,ae) - f[4+(w-1)*4, i])/4)
            CM_Dyn = CM_st + CL_Dyn*(ast_1 - ast_2) - (pi/2)*Tu*alphadot(tvec[i])


            
            coefficient_matrix[1+(w-1)*4 , i] = alpha(tvec[i]) #The coefficients and angle of attack are pushed into the matric made at the beginning of the function
            coefficient_matrix[2+(w-1)*4 , i] = CL_Dyn
            coefficient_matrix[3+(w-1)*4 , i] = CD_Dyn
            coefficient_matrix[4+(w-1)*4 , i] = CM_Dyn
            
        end
    end

    return coefficient_matrix
end