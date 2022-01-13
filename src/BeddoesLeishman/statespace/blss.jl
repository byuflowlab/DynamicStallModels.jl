
function ffun(alpha, alpha1, S)
    if alpha<=alpha1
        return 1.0-0.3*exp((alpha-alpha1)/S[1])
    else
        return 0.04 + 0.66*exp((alpha1-alpha)/S[2]) #One paper has a plus here, another has a minus. 
    end
end

function prepenvironment(; c=0.1, M=0.3, a=343.3, amp=10.0, shift=10.0, k=0.1)
    v = M*a
    omega = k*2*v/c
    
    function U(t)
        return M*a
    end

    function Udot(t)
        return 0
    end

    function alpha(t)
        alf = shift + amp*sin(omega*t)
        return alf*(pi/180)
    end

    function alphadot(t)  
        alfd = amp*omega*cos(omega*t)
        return alfd*(pi/180)
    end

    function alphaddot(t)
        alfdd = -amp*(omega^2)*sin(omega*t)
        return alfdd*(pi/180)
    end
    
    return U, Udot, alpha, alphadot, alphaddot
end


function Heavi(t)
    if t>0
        return 1
    else
        return 0
    end
end

function states!(dx, x, p, t)
    ### Unpack p
    V, Vdot, alpha, alphadot, c, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv = p

    
    ### Environmental values
    M = V(t)/a #Mach number
    beta = sqrt(1-(M^2)) #Compressibility factor #TODO: This will need to be updated to handle if someone hands it a mach number greater than 1. 
    q(t) = alphadot(t)*c/V(t) #Pitch rate

    ### Constants
#     conv = (c/(2*V(t)))
    tau_p = T_p*conv(t)
    tau_f = T_f*conv(t)
    tau_vl = T_vl*conv(t) #Todo: Should I multiply this by conv
    tau_v = T_v*conv(t)

    # println(tau_vl)

    tau_l = c/a #Todo. I need to find out if this should be the definition just after equation 4 in leishman 1990, or on the top right page of 839 of the same paper. This is defined in the nomenclature section. 
    k_alpha = 1/((1-M) + (pi*beta*(M^2)*((A[1]*b[1])+(A[2]*b[2]))))
    k_q = 1/((1-M) + (2*pi*beta*(M^2)*((A[1]*b[1])+(A[2]*b[2]))))

    ### Update the attached flow states (The first four)
    dx[1] = (2*V(t)/c)*(beta^2)*(-b[1]*x[1]) + alpha(t) + q(t)/2 #High frequency shed vortex
    dx[2] = (2*V(t)/c)*(beta^2)*(-b[2]*x[2]) + alpha(t) + q(t)/2 #Low frequency shed vortex
    dx[3] = (-1/(k_alpha*tau_l))*x[3] + alpha(t) #Impulse due to angle of attack
    dx[4] = (-1/(k_q*tau_l))*x[4] + q(t) #Impulse due to pitch rate. 

    ### Calculate the attached lift
    Ccn = dCndalpha*(2*V(t)/c)*(beta^2)*((A[1]*b[1]*x[1])+(A[2]*b[2]*x[2])) #Circulatory lift
    Cina = (4/M)*dx[3] #Normal force due to angle of attack #Todo: Should I calculate this before I update dx? 
    Cinq = (1/M)*dx[4] #Normal force due to pitch rate #Todo: Same as above. I can look at the indicial approach to get an idea. 
    Cpn = Ccn + Cina + Cinq #Attached flow lift
    # println("Cpn: ", Cpn)

    ### Calculate the dynamic flow states (the second four)
    dx[5] = (-1/tau_p)*x[5] + (Cpn/tau_p) #Delayed attached lift.
    # println("tau_p: ", tau_p)

    alpha_f = x[5]/dCndalpha #Equivalent angle of attack
    fp = ffun2(alpha_f, alpha1, S.*conv(t)) #Equivalent separation point
    # println("fp: ", fp)
    # if fp<0
    #     println("")
    #     println("alpha_f: ", alpha_f)
    #     println("fp: ", fp)
    # end
    dx[6] = (-1/tau_f)*x[6] + (fp/tau_f) #Delayed equivalent separation point

    dx[7] = (V(t)/(3*c))*Heavi(x[5]-Cn1)*Heavi(2*tau_vl - x[7]) - V(t)*Heavi(Cn1-x[5])*Heavi(x[7])/c #Position of the separation point
    
    Ccndot = dCndalpha*(2*Vdot(t)/c)*(1-(M^2))*((A[1]*b[1]*x[1])+(A[2]*b[2]*x[2])) + dCndalpha*(2*V(t)/c)*(1-(M^2))*((A[1]*b[1]*dx[1])+(A[2]*b[2]*dx[2])) - dCndalpha*(4*(M^2)/c)*Vdot(t)*((A[1]*b[1]*x[1])+(A[2]*b[2]*x[2])) #Derivative of the circulatory normal force with respect to time. 
    Cvdot = (Ccndot*(1-((1+sqrt(x[6]))^2)/4) - Ccn*(1 + sqrt(x[6]))*dx[6]/(4*sqrt(x[6])))*Heavi(2*tau_vl - x[7])*Heavi(x[5]-Cn1) #Derivative of the vortex lift contribution with respect to time.
    
#     if Cvdot>0.001
#         println(Cvdot)
#     end
    # println("") 
    # println("Ccn: ", Ccn)
    # println("Ccndot: ", Ccndot)
    # println("Cvdot: ", Cvdot)
    dx[8] = (-1/tau_v)*x[8] + (Cvdot/tau_v) #Vortex lift contribution
end


function parsesolution(sol, p; eta=0.95)
    function extractdata(sol)
        t = sol.t
        u = reduce(hcat, sol.u)'
        n,m = size(u)
        du = zeros(n,m)

        for j=1:m
            du[:,j] = gradient(t, u[:,j], t)
        end
        return t, u, du
    end

    V, Vdot, alpha, alphadot, c, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv = p

    t, x, dx = extractdata(sol)

    n = length(t)

    #Initialize arrays
    Ccn = zeros(n)
    Cina = zeros(n)
    Cinq = zeros(n)
    Cpn = zeros(n)

    alphaE = zeros(n)
    Cnf = zeros(n)
    Cnt = zeros(n)
    
    Cfc = zeros(n)
    Cdn = zeros(n)

    for i = 1:n
        ### Calculate constants
        M = V(t[i])/a
        beta = sqrt(1-(M^2))

        ### Calculate attached lift
        Ccn[i] = dCndalpha*(2*V(t[i])/c)*(beta^2)*((A[1]*b[1]*x[i,1])+(A[2]*b[2]*x[i,2])) #Circulatory lift
        Cina[i] = (4/M)*dx[i,3] #Normal force due to angle of attack  
        Cinq[i] = (1/M)*dx[i,4] #Normal force due to pitch rate 

        Cpn[i] = Ccn[i] + Cina[i] + Cinq[i] #Attached normal force

        alphaE[i] = (beta^2)*(2*V(t[i])/c)*((A[1]*b[1]*x[i,1])+(A[2]*b[2]*x[i,2])) #Effective angle of attack

        Cnf[i] = dCndalpha*(((1+sqrt(x[i,6]))^2)/4)*alphaE[i] + Cina[i] + Cinq[i] #Nonlinear dynamic normal force

        Cnt[i] = Cnf[i] + x[i,8] #Total dynamic normal force
        
        Cfc[i] = eta*dCndalpha*(alphaE[i]^2)*sqrt(x[i,6])
        Cdn[i] = Cnt[i]*sin(alpha(t[i])) - Cfc[i]*cos(alpha(t[i]))
    end
    return t, x, dx, Cnt, Cpn, Cdn, Cfc
end