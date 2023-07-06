export Riso

struct Riso<: DSModel
    detype::DEType
    A
    b
    T
end


function numberofstates(dsmodel::Riso)
    return 4
end

function numberofparams(dsmodel::Riso)
    return 6
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

function initialize(dsmodel::Riso, airfoil::Airfoil, tvec, y)
    if isa(dsmodel.detype, Indicial)
        states = [0.0, 0.0, 0.0, 0.0]

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


function update_states!(model::Riso, airfoil::Airfoil, oldstate, newstate, y, dt)
    U_new, Udot_new, alpha_new, alphadot_new = y
    #find out a way to get the previous parameter values
    x1_old , x2_old, x3_old, x4_old = oldstate


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


    P1 = b1*(U_new + U_old)/c + (Udot_new + Udot_old)/(U_new + U_old)
    P2 = b2*(U_new + U_old)/c + (Udot_new + Udot_old)/(U_new + U_old)
    P3 = 1/Tp
    P4 = 1/Tf

    Q1 = (U_old*alpha_old + U_new*alpha_new)*(b1*A1)/(c)
    Q1 = (U_old*alpha_old + U_new*alpha_new)*(b2*A2)/(c)
    Q3 = (0.5/Tp)*(dcldalpha*(ae_old + ae_new - 2*alpha0) + (pi*c*0.5)*((alphadot_old/U_old + alphadot_new/U_new)))

    x1_new = exp(-1*P1*dt)*x1_old + (Q1/P1)*(1 - exp(-1*P1*dt))
    x2_new = exp(-1*P2*dt)*x2_old + (Q2/P2)*(1 - exp(-1*P2*dt))
    x3_new = exp(-1*P3*dt)*x3_old + (Q3/P3)*(1 - exp(-1*P3*dt))

    Q4 = (0.5/Tf)*(separationpoint(airfoil, x3_old/(dcldalpha+alpha0)) + separationpoint(airfoil, x3_new/(dcldalpha+alpha0)))

    x4_new = exp(-1*P4*dt)*x4_old + (Q4/P4)*(1 - exp(-1*P4*dt))

    newstate[1] = x1_new
    newstate[2] = x2_new
    newstate[3] = x3_new
    newstate[4] = x4_new
end


function get_loads(model::Riso, airfoil::Airfoil, states, y)
    x1, x2, x3, x4 = states

    U, _, alpha, alphadot = y

    alpha0 = airfoil.alpha0
    c = airfoil.c
    dcldalpha = airfoil.dcldalpha
    CD0 = airfoil.cd(alpha0)
    CM0 = airfoil.cm(alpha0)

    A1 = model.A[1]
    A2 = model.A[2]


    ae = alpha*(1 - A1 - A2) + x1 + x2
    Tu = c/(U*2)
    CL_fs = hansen_fully_sep(airfoil , ae)
    CD_st = airfoil.cd(ae)
    CM_st = airfoil.cm(ae)
    ast_1 = (airfoil.cm(x4 - CM0)/airfoil.cl(x4) 
    ast_2 = (airfoil.cm(separationpoint(airfoil, ae)) - CM0)/airfoil.cl(separationpoint(airfoil, ae)) 

    CL_Dyn = dcldalpha*(ae - alpha0)*x4 + CL_fs*(1 - x4) + pi*Tu*alphadot
    CD_Dyn = CD_st + (alpha - ae)*CL_Dyn + (CD_st - CD0)*((sqrt(separationpoint(airfoil, ae)) - sqrt(x4))/2 - (separationpoint(airfoil,ae) - x4)/4)
    CM_Dyn = CM_st + CL_Dyn*(ast_1 - ast_2) - (pi/2)*Tu*alphadot

    return CL_Dyn, CD_Dyn, CM_Dyn

end

function get_loads!(model::Riso, airfoil::Airfoil, states, loads, y)
    loads .= get_loads(model::Riso, airfoil::Airfoil, states, y)
end


function state_rates!(model::Riso , airfoil::Airfoil , dx, x, y, t)
    
    U, Udot, alpha, alphadot = evaluate_environment(y,t)


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

    ae = alpha*(1-A1-A2) + x[1] + x[2]

    fst_val = separationpoint(airfoil, (alpha0 + x[3]/dcldalpha))

    dx[1] = b1*A1*alpha/Tu - x[1]*(b1+c*Udot/(2*U^2))/Tu
    dx[2] = b2*A2*alpha/Tu - x[2]*(b2 + c*Udot/(2*U^2))/Tu
    dx[3] = (dcldalpha*(ae-alpha0) + pi*Tu*alphadot)/Tp - x[3]/Tp
    dx[4] = fst_val/Tf - x[4]/Tf
end


function hansen_fully_sep(airfoil , alpha)
    CL_st = airfoil.cl(alpha)
    dcldalpha = airfoil.dcldalpha
    alpha0 = airfoil.alpha0
    delta_alpha = (alpha-alpha0)

    fst = (2*sqrt(CL_st/(dcldalpha*(delta_alpha))) - 1)^2

    CL_fs = (CL_st - dcldalpha*delta_alpha*fst)/(1-fst)

    return CL_fs
end

export parsesolution

function parsesolution(model::Riso, airfoils::AbstractVector{<:Airfoil}, sol::ODESolution, y)
    f = Array(sol)
    tvec = sol.t


    coefficient_matrix = zeros(4*length(airfoils) , length(tvec))
    fs_matrix = zeros(1,length(tvec))
    Cl_matrix = zeros(1,length(tvec))
    sp_matrix = zeros(1,length(tvec))

    for w in 1:length(airfoils)

        airfoil = airfoils[w]

        alpha0 = airfoil.alpha0
        c = airfoil.c
        U = y[1+(w-1)*4]
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



            CL_Dyn = dcldalpha*(ae - alpha0)*f[4+(w-1)*4, i] + CL_fs*(1 - f[4+(w-1)*4, i]) + pi*Tu*alphadot(tvec[i])
            CD_Dyn = CD_st + (alpha(tvec[i]) - ae)*CL_Dyn + (CD_st - CD0)*((sqrt(separationpoint(airfoil, ae)) - sqrt(f[4+(w-1)*4, i]))/2 - (separationpoint(airfoil,ae) - f[4+(w-1)*4, i])/4)
            CM_Dyn = CM_st + CL_Dyn*(ast_1 - ast_2) - (pi/2)*Tu*alphadot(tvec[i])


            
            coefficient_matrix[1+(w-1)*4 , i] = alpha(tvec[i])
            coefficient_matrix[2+(w-1)*4 , i] = CL_Dyn
            coefficient_matrix[3+(w-1)*4 , i] = CD_Dyn
            coefficient_matrix[4+(w-1)*4 , i] = CM_Dyn
            fs_matrix[1,i] = hansen_fully_sep(airfoil, alpha(tvec[i]))
            Cl_matrix[1,i] = airfoil.cl(alpha(tvec[i]))
            sp_matrix[1,i] = separationpoint(airfoil, alpha(tvec[i]))
        end
    end

    return coefficient_matrix, fs_matrix, Cl_matrix, sp_matrix
end