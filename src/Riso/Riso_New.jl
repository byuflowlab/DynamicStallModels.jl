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