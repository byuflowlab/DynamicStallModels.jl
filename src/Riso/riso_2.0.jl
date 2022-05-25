
function seperationpoint(alpha, afm, afp, clfit, dcldalpha, alpha0)
    if alpha<afm #The theory says that if the angle of attack is less (or greater than) than the seperation point aoa, then the function should return 0. 
        return typeof(alpha)(0)
    elseif alpha>afp
        return typeof(alpha)(0)
    end
    
    cl_static = clfit(alpha)
    cl_linear = dcldalpha*(alpha-alpha0)
    f = (2*sqrt(cl_static/(cl_linear))-1)^2

    if f>1
        return typeof(alpha)(1)
    elseif isnan(f)
        return typeof(alpha)(1)
    end
    return f
end



function riso_state_rates!(dx, x, U, Udot, alpha, alphadot, c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp)

    Tu = c/(2*U)
    alpha_E = alpha*(1 - A1 - A2) + x[1] + x[2]
    alpha_f = x[3]/dcldalpha + alpha0
    fst = seperationpoint(alpha_f, afm, afp, clfit, dcldalpha, alpha0)

    dx[1] = -x[1]*(b1 + c*Udot/(2*(U^2)))/Tu + b1*A1*alpha/Tu
    dx[2] = -x[2]*(b2 + c*Udot/(2*(U^2)))/Tu + b2*A2*alpha/Tu
    dx[3] = -x[3]/Tp + (dcldalpha*(alpha_E - alpha0) + pi*Tu*alphadot)/Tf
    dx[4] = -x[4]/Tf + fst/Tf
end

function riso_ode!(dx, x, p, t)
    U, Udot, alpha, alphadot, c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp = p

    riso_state_rates!(dx, x, U(t), Udot(t), alpha(t), alphadot(t), c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp)
end

function Clfs(alpha, clfit, dcldalpha, alpha0, afm, afp)
    fst = seperationpoint(alpha, afm, afp, clfit, dcldalpha, alpha0)

    
    clst = clfit(alpha)
    # @show fst
    if fst>=1 #Todo: This is supposed to happen automagically. 
        return clst/2
    else
        return (clst - dcldalpha*(alpha - alpha0)*fst)/(1-fst)
    end
    
end

function riso_coefficients(x, U, alpha, alphadot, c, clfit, cdfit, dcldalpha, Cd0, alpha0, A1, A2, afm, afp)

    Tu = c/(2*U)
    alpha_E = alpha*(1 - A1 - A2) + x[1] + x[2]

    # alpha_f = x[3]/dcldalpha + alpha0
    # fst = seperationpoint(alpha_f, afm, afp, clfit, dcldalpha, alpha0)

    ### Calculate Clfs #Todo: Should Clfs use fst calculated from alpha_f or alpha? I'm going to assume it is based off of alpha, but recognize that this could be wrong. 
    cl_fs = Clfs(alpha_E, clfit, dcldalpha, alpha0, afm, afp)
    # @show cl_fs

    Cl_dyn = dcldalpha*(alpha_E - alpha0)*x[4] + cl_fs*(1-x[4]) + pi*Tu*alphadot

    cdst = cdfit(alpha_E)
    fst = seperationpoint(alpha_E, afm, afp, clfit, dcldalpha, alpha0)
    t1 = (sqrt(fst) - sqrt(x[4]))/2 #TODO: Here is the square root that causes problems when the fourth state is negative. 
    t2 = (fst - x[4])/4
    Cd_dyn = cdst + (alpha - alpha_E)*Cl_dyn + (cdst - Cd0)*(t1 - t2)
    return Cl_dyn, Cd_dyn
end

function parsesolution(sol, p, cdfit, Cd0)

    ### Unpack
    x = Array(sol)'
    t = sol.t
    U, Udot, alpha, alphadot, c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp = p

    nt = length(t)

    clvec = zeros(nt)
    cdvec = zeros(nt)

    ### run through the time steps and calculate the dynamic lift and drag (based on the states)
    for i = 1:nt
        ti = t[i]
        clvec[i], cdvec[i] = riso_coefficients(x[i,:], U(ti), alpha(ti), alphadot(ti), c, clfit, cdfit, dcldalpha, Cd0, alpha0, A1, A2, afm, afp)
    end

    return clvec, cdvec, t
end