



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










function states_liftonly!(du, u, p, t)
    ### Unpack p
    V, alpha, alphadot, c, dcldalpha, alpha1, alpha0, Cn1, A, b, a, Tp, Tf, Tv = p

    q(t) = alphadot(t)*c/V(t)

    ### Calculate constants
    M = V(t)/a #Mach number
    beta = sqrt(1-(M^2)) #Compressibility factor

    ### Calculate time constants
    Tl = c/a
    k_alpha = 1/((1-M) + (pi*beta*(M^2)*(A[1]*b[1] + A[2]*b[2]))) #I wonder if any of these should include slope. -> Nope, they are included in the normal force equations. 

    k_q = 1/((1-M) + (2*pi*beta*(M^2)*(A[1]*b[1] + A[2]*b[2])))

    ### Calculate matrix coefficients
    a11 = -2*b[1]*V(t)*(beta^2)/c
    a22 = -2*b[2]*V(t)*(beta^2)/c
    a32 = -1/(k_alpha*Tl)
    a44 = -1/(k_q*Tl)
    
    ### Calculate the attached state rates 
    #Shed wake terms
    du[1] = u[1]*a11 + alpha(t) + q(t)/2
    du[2] = u[2]*a22 + alpha(t) + q(t)/2

    #noncirculatory angle of attack contribution
    du[3] = u[3]*a32 + alpha(t)

    #noncirculatory pitch rate contribution
    du[4] = u[4]*a44 + q(t)

    ### Calculate the attached valules
    #Circulatory normal force
    Ccn = dcldalpha*(2*V(t)/c)*(beta^2)*(A[1]*b[1]*u[1] + A[2]*b[2]*u[2])

    #noncirculatory normal force due to angle of attack
    Cina = 4*du[3]/M

    #noncirculatory normal force due to pitch rate
    Cinq = du[4]/M

    #Total attached flow normal force
    Cpn = Ccn + Cina + Cinq 

    ### Calculate the dynamic state rates
    #Note: These state rates are labeled differently in the notes. It's because I'm omitting the moment states

    #delayed attached lift
    du[5] = (Cpn-u[5])/Tp #-u[5]/Tp + Cpn/Tp #sv9
    
    # Delayed effective angle of attack
    alpha_f = u[5]/dcldalpha

    # Unsteady TE separation point, f''
    du[6] = -u[6]/Tf + separationpoint_BL(alpha_f, alpha1, alpha0)/Tf #sv10 

    tau_v = 2*Tv

    # Separation point location time constant, τ
    du[7] = V(t)*Heavi(u[5]-Cn1)*Heavi(tau_v-u[7])/(3*c) - V(t)*Heavi(Cn1-u[5])*Heavi(u[7])/c #Note: This is my own creation. Especially the stuff after the negative sign. I figure it will take some time to reattach fully, but I do know that there can be multiple vortices attached... so it probably should just reset to 0. Or maybe reset faster.  

    # Change in lift due to vortex position
    Ccndot = dcldalpha*(2*V(t)/c)*(beta^2)*(A[1]*b[1]*du[1] + A[2]*b[2]*du[2]) 
    Cvdot = (Ccndot*(1- ((1+sqrt(u[6]))^2)/4) - Ccn*(1+sqrt(u[6]))*du[6]/(4*sqrt(u[6])))*Heavi(tau_v-u[7])*Heavi(u[5]-Cn1) #Note: If u[6]<=0, then we be screwed. (u[6]=f'') 

    # Vortex lift contribution
    du[8] = -u[8]/Tv + Cvdot/Tv  
end

function states!(dx, x, p, t)
    ### Unpack p
    V, Vdot, alpha, alphadot, c, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv = p

    
    ### Environmental values
    M = V(t)/a #Mach number
    beta = sqrt(1-(M^2)) #Compressibility factor #TODO: This will need to be updated to handle if someone hands it a mach number greater than 1. 
    q(t) = alphadot(t)*c/V(t) #Pitch rate

    ### Constants
    #conv = (c/(2*V(t)))
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
    fp = separationpoint_BL(alpha_f, alpha1, S.*conv(t)) #Equivalent separation point
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

function states_attachedflow!(du, u, p, t)
    #Leishman 1990
    ### Unpack p

    ### Calculate constants
    V = sqrt(u(t)^2 + v(t)^2) #TODO: I probably won't need v, but I'll leave it for now. I just need someway to make sure that I can fully  account for the inflow speed and angle, and the commanded pitch of the airfoil. It might be a good idea to have a function that converts from u, v , and pitch to U and alpha. 
    M = V/a
    beta = sqrt(1-(M^2))

    ### Calculate time constants
    Tl = c/a
    k_alpha = 1/((1-M) + (pi*beta*(M^2)*(A[1]*b[1] + A[2]*b[2]))) #I wonder if any of these should include slope. -> Nope, they are included in the normal force equations. 

    k_q = 1/((1-M) + (2*pi*beta*(M^2)*(A[1]*b[1] + A[2]*b[2])))

    k_alpha_m = (A[3]*b[3] + A[4]*b[4])/(b[3]*b[4]*(1-M))

    k_q_m = 7/(15*(1-M) + (3*pi*beta*(M^2)*b[5]))

    ### Calculate matrix coefficients
    a11 = -2*b[1]*V*(beta^2)/c
    a22 = -2*b[2]*V*(beta^2)/c
    a32 = -1/(k_alpha*Tl)
    a44 = -1/(k_q*Tl)
    a55 = -1/(b[3]*k_alpha_m*Tl)
    a66 = -1/(b[4]*k_alpha_m*Tl)
    a77 = -2*b[5]*V*(beta^2)/c
    a88 = -1/(k_q_m*Tl)

    ### Calculate the state rates 
    # I decided to leave the state rates as a 1D array rather than a matrix. A matrix would have been easier to code... but I'm not sure if it would have been more computationally efficient. 
    du[1] = u[1]*a11 + alpha(t) + q(t)/2
    du[2] = u[2]*a22 + alpha(t) + q(t)/2
    du[3] = u[3]*a32 + alpha(t)
    du[4] = u[4]*a44 + q(t)
    du[5] = u[5]*a55 + alpha(t)
    du[6] = u[6]*a66 + alpha(t)
    du[7] = u[7]*a77 + q(t)
    du[8] = u[8]*a88 + q(t)

end

function fullstates!(dx, x, p, t) #Todo: All of these equations should be checked. 
    ### Unpack p
    V, Vdot, alpha, alphadot, c, dCndalpha, alpha1, Cn1, Cmo, A, b, S, K, m, T_p, T_f, T_vl, T_v, a = p #Note: I do use Vdot. 

    ### Environmental values
    M = V(t)/a #Mach number
    beta = sqrt(1-(M^2)) #Compressibility factor #TODO: This will need to be updated to handle if someone hands it a mach number greater than 1. 
    q(t) = alphadot(t)*c/V(t) #Pitch rate

    ### Constants
    L = 2*V(t)/c
    tau_p = T_p/L
    tau_f = T_f/L
    tau_vl = T_vl/L
    tau_v = T_v/L

    # println(tau_vl)

    tau_l = c/a #Todo. I need to find out if this should be the definition just after equation 4 in leishman 1990, or on the top right page of 839 of the same paper. This is defined in the nomenclature section. 
    k_alpha = 1/((1-M) + (pi*beta*(M^2)*((A[1]*b[1])+(A[2]*b[2]))))
    k_q = 1/((1-M) + (2*pi*beta*(M^2)*((A[1]*b[1])+(A[2]*b[2]))))
    k_alpha_m = ((A[3]*b[4])+(A[4]*b[3]))/(b[3]*b[4]*(1-M))
    k_q_m = 7/((15*(1-M))+(3*pi*beta*(M^2)*b[5]))

    ### Update the attached flow states (The first four)
    dx[1] = L*(beta^2)*(-b[1]*x[1]) + alpha(t) + q(t)/2 #High frequency shed vortex
    dx[2] = L*(beta^2)*(-b[2]*x[2]) + alpha(t) + q(t)/2 #Low frequency shed vortex
    dx[3] = (-1/(k_alpha*tau_l))*x[3] + alpha(t) #Impulse due to angle of attack
    dx[4] = (-1/(k_q*tau_l))*x[4] + q(t) #Impulse due to pitch rate. 
    dx[5] = (-1/(b[3]*k_alpha_m*tau_l))*x[5] + alpha(t)
    dx[6] = (-1/(b[4]*k_alpha_m*tau_l))*x[6] + alpha(t)
    dx[7] = -b[5]*(beta^2)*L*x[7] + q(t)
    dx[8] = (-1/(k_q_m*tau_l))*x[8] + q(t)

    ### Calculate the attached lift
    Ccn = dCndalpha*L*(beta^2)*((A[1]*b[1]*x[1])+(A[2]*b[2]*x[2])) #Circulatory lift
    Cina = (4/M)*dx[3] #Normal force due to angle of attack #Todo: Should I calculate this before I update dx? 
    Cinq = (1/M)*dx[4] #Normal force due to pitch rate #Todo: Same as above. I can look at the indicial approach to get an idea. 
    Cpn = Ccn + Cina + Cinq #Attached flow lift
    # println("Cpn: ", Cpn)


    ### Calculate the dynamic flow states (the second four)
    dx[9] = (-1/tau_p)*x[9] + (Cpn/tau_p) #Delayed attached lift.
    # println("tau_p: ", tau_p)

    alpha_f = x[9]/dCndalpha #Equivalent angle of attack
    fp = separationpoint_BL(alpha_f, alpha1, S./L) #Equivalent separation point #S.*conv(t)
    # println("fp: ", fp)
    # if fp<0
    #     println("")
    #     println("alpha_f: ", alpha_f)
    #     println("fp: ", fp)
    # end
    dx[10] = (-1/tau_f)*x[10] + (fp/tau_f) #Delayed equivalent separation point

    # dx[11] = (V(t)/(2*c))*Heavi(x[9]-Cn1)*Heavi(2*tau_vl - x[11]) - V(t)*Heavi(Cn1-x[5])*Heavi(x[7])/c #Position of the separation point
    # dx[11] = (V(t)/(2*c))*Heavi(x[9]-Cn1)*Heavi(2*tau_vl - x[11]) #I want to try without my diminishing term
    # dx[11] = (V(t)/(3*c))*Heavi(x[9]-Cn1) #I'm not sure that this should cut off when it passes the 2tau value, what if I want a second vortex?
    dx[11] = (2*V(t)/c)*Heavi(x[9]-Cn1) #The dymore has delta t * v/b... which would imply that the 2 is on top, not bottom. 
    # dx[11] = Heavi(x[9]-Cn1)/3 #Dimensional vortex location (seconds) I'm not sure if it makes sense to do this. 
    if x[9]<Cn1 #TODO: This can easily be rewritten as x[11] = x[11]*Heavi(abs(x[9])-Cn1) #Which means that tau will be reset if Cn' is less then Cn1
        x[11]=0
    end
    
    Ccndot = dCndalpha*(2*Vdot(t)/c)*(1-(M^2))*((A[1]*b[1]*x[1])+(A[2]*b[2]*x[2])) + dCndalpha*(2*V(t)/c)*(1-(M^2))*((A[1]*b[1]*dx[1])+(A[2]*b[2]*dx[2])) - dCndalpha*(4*(M^2)/c)*Vdot(t)*((A[1]*b[1]*x[1])+(A[2]*b[2]*x[2])) #Derivative of the circulatory normal force with respect to time. 

    # Cvdot = (Ccndot*(1-((1+sqrt(x[10]))^2)/4) - Ccn*(1 + sqrt(x[10]))*dx[10]/(4*sqrt(x[10])))*Heavi(2*tau_vl - x[11])*Heavi(x[9]-Cn1) #*(c/a) #*(2*V(t)/c) #Derivative of the vortex lift contribution with respect to time. #TODO: There is something going on here. I don't think the units quite work out. I tried correcting the units by multiplying by 2V/c, but that gave me this huge jump. And I don't think there is supposed to be this huge jump. At least, that's what the paper's solutions show. Thus there is something up. I also tried multiplying by c/a... which doesn't help because it makes the forcing function even smaller. 
    # Cvdot = Ccn*(1-(((1+sqrt(x[10]))/2)^2))*Heavi(2*tau_vl - x[11])*Heavi(x[9]-Cn1) #Note: Cvn dropped even more. 
    # Cvdot = (Ccndot*(1-((1+sqrt(x[10]))^2)/4) - Ccn*(1 + sqrt(x[10]))*dx[10]/(4*sqrt(x[10])))*Heavi(2*tau_vl - x[11]) - Ccn*(1-((1+sqrt(x[10]))^2)/4)*dx[11]*diracdelta(2*tau_vl-x[11]) #I'm comparing a united time constant against a unitless state.
    # Cvdot = (Ccndot*(1-((1+sqrt(x[10]))^2)/4) - Ccn*(1 + sqrt(x[10]))*dx[10]/(4*sqrt(x[10])))*Heavi(2*T_vl - x[11]) - Ccn*(1-((1+sqrt(x[10]))^2)/4)*dx[11]*diracdelta(2*T_vl-x[11]) 
    #Comparing against unitless time constants made a huge difference, now the results are in the right direction, but much too huge. I wonder if it's that the rate of tau is too large? 
    # Cvdot = (Ccndot*(1-((1+sqrt(x[10]))^2)/4) - Ccn*(1 + sqrt(x[10]))*dx[10]/(4*sqrt(x[10])))*Heavi(T_vl - x[11]) - Ccn*(1-((1+sqrt(x[10]))^2)/4)*dx[11]*diracdelta(T_vl-x[11]) #That didn't make too big of a change. 
    Cvdot = (Ccndot*(1-((1+sqrt(x[10]))^2)/4) - Ccn*(1 + sqrt(x[10]))*dx[10]/(4*sqrt(x[10])))*Heavi(2*T_vl - x[11]) #Note: This was the working one. 
    # Cvdot = (Ccndot*(1-((1+sqrt(x[10]))^2)/4) - Ccn*(1 + sqrt(x[10]))*dx[10]/(4*sqrt(x[10])))*Heavi(2*tau_vl - x[11])
    # Cvdot = (Ccndot*(1-((1+sqrt(x[10]))^2)/4) - Ccn*(1 + sqrt(x[10]))*dx[10]/(4*sqrt(x[10])))*Heavi(T_vl - x[11]) - Ccn*(1-((1+sqrt(x[10]))^2)/4)*dx[11]*dirac(T_vl-x[11])
    # Cvdot = (Ccndot*(1-((1+sqrt(x[10]))^2)/4) - Ccn*(1 + sqrt(x[10]))*dx[10]/(4*sqrt(x[10])))*Heavi(x[11] - 2*T_vl)
    


    #     if Cvdot>0.001
    #         println(Cvdot)
    #     end
    # println("") 
    # println("Ccn: ", Ccn)
    # println("Ccndot: ", Ccndot)
    # println("Cvdot: ", Cvdot)
    # dx[12] = (-1/tau_v)*x[12] + Cvdot/(tau_v) #Vortex lift contribution
    dx[12] = (-1/tau_v)*x[12] + Cvdot/(T_v) #Vortex lift contribution
    # dx[12] = (-1/tau_v)*x[12] + Cvdot/(tau_v*(2*V(t)/c)) #Trying by using a lame way to derive Cv w.r.t nondimensional time. 
    dx[13] = (-2/tau_f)*x[13] + 2*fp/tau_f

end

function parsesolution_liftonly(sol, p; eta=0.95)
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

    q(t) = alphadot(t)*c/V(t)

    V, alpha, alphadot, c, dcldalpha, alpha1, alpha0, Cn1, A, b, a, Tp, Tf, Tv = p

    tvec, u, du = extractdata(sol)
    n = length(tvec)

    #Initialize output vectors. 
    Cl = zeros(n)
    Cd = zeros(n)

    Cn = zeros(n)
    Ccn = zeros(n)
    Cina = zeros(n)
    Cinq = zeros(n)
    Cpn = zeros(n)
    Cpn2 = zeros(n)

    alphaE = zeros(n)
    Cf = zeros(n)
    Cdn = zeros(n)
    Ccndot = zeros(n)
    K = zeros(n)

    Cfc = zeros(n)
    Cdc = zeros(n)
    Cdd = zeros(n)
    

    for i=1:n
        ti = tvec[i]

        ### Calculate constants
        M = V(ti)/a #Mach number
        beta = sqrt(1-(M^2)) #Compressibility factor

        ### Calculate time constants
        Tl = c/a
        k_alpha = 1/((1-M) + (pi*beta*(M^2)*(A[1]*b[1] + A[2]*b[2])))  
        k_q = 1/((1-M) + (2*pi*beta*(M^2)*(A[1]*b[1] + A[2]*b[2])))

        ### Attached flow
        #Circulatory normal force
        Ccn[i] = dcldalpha*(2*V(ti)/c)*(beta^2)*(A[1]*b[1]*u[i,1] + A[2]*b[2]*u[i,2])

        #noncirculatory normal force due to angle of attack
        Cina[i] = 4*du[i,3]/M 

        #noncirculatory normal force due to pitch rate
        Cinq[i] = du[i,4]/M 

        #Total attached flow normal force
        Cpn[i] = Ccn[i] + Cina[i] + Cinq[i]

        ### Dynamic Stall  
        alphaE[i] = (beta^2)*(2*V(ti)/c)*(A[1]*b[1]*u[i,1] + A[2]*b[2]*u[i,2]) #Equivalent angle of attack (due to shed wake terms)
        K[i] = ((1+sqrt(u[i,6]))^2)/4
        Cf[i] = dcldalpha*K[i]*alphaE[i] + Cina[i] + Cinq[i] 
        Cdn[i] = Cf[i] + u[i,8] #Total normal force - dynamic stall conditions

        ### Drag
        Cfc[i]  = eta*dcldalpha*((alphaE[i])^2)*sqrt(u[i,6]) #Nonlinear component

        Cdc[i] = Cfc[i] #It think there should be a vortex component.


        Cdd[i] = Cdn[i]*sin(alpha(ti)) - Cdc[i]*cos(alpha(ti))

    end
    return Cdn, Cdd, Cpn, u, du
end

##### Parsing the solution while trying several different calculations for drag. -> Which would suggest that the lift portion of the DS model was working quite well. 

# function parsesolution(sol, p; eta=0.95)
#     function extractdata(sol)
#         t = sol.t
#         u = reduce(hcat, sol.u)'
#         n,m = size(u)
#         du = zeros(n,m)

#         for j=1:m
#             du[:,j] = gradient(t, u[:,j], t)
#         end
#         return t, u, du
#     end

#     q(t) = alphadot(t)*c/V(t)

#     V, alpha, alphadot, c, dcldalpha, alpha1, alpha0, Cn1, A, b, a, Tp, Tf, Tv = p

#     tvec, u, du = extractdata(sol)
#     n = length(tvec)

#     Cl = zeros(n)
#     Cd = zeros(n)

#     Cn = zeros(n)
#     Ccn = zeros(n)
#     Cina = zeros(n)
#     Cinq = zeros(n)
#     Cpn = zeros(n)
#     Cpn2 = zeros(n)

#     alphaE = zeros(n)
#     Cf = zeros(n)
#     Cdn = zeros(n)
#     Ccndot = zeros(n)
#     K = zeros(n)

#     Cc = zeros(n)
#     Cdc = zeros(n)

#     Cdd = zeros(n)
#     Cdd2 = zeros(n)
#     Cdd3 = zeros(n)
#     Cdd4 = zeros(n)
#     Cdd5 = zeros(n)
#     Cdd6 = zeros(n)
#     Cdd7 = zeros(n)
#     Cdd8 = zeros(n)
#     Cdd9 = zeros(n)
#     Cdd10 = zeros(n)
#     Cdd11 = zeros(n)
#     Cdd12 = zeros(n)
#     Cdd13 = zeros(n)
#     Cdd14 = zeros(n)
#     Cdd15 = zeros(n)
#     Cdd16 = zeros(n)
#     Cdd17 = zeros(n)
#     Cdd18 = zeros(n)
#     Cdd19 = zeros(n)

#     Cfc = zeros(n)
#     Cfc2 = zeros(n)
#     Cfc3 = zeros(n)
#     Cfc4 = zeros(n)
#     Cfc5 = zeros(n)
#     Cfc6 = zeros(n)
#     Cfc7 = zeros(n)
#     Cfc8 = zeros(n)
#     Cfc9 = zeros(n)
#     Cfc10 = zeros(n)
#     Cfc11 = zeros(n)
#     Cfc12 = zeros(n)
#     Cfc13 = zeros(n)
#     Cfc14 = zeros(n)
#     Cfc15 = zeros(n)
#     Cfc16 = zeros(n)
#     Cfc17 = zeros(n)
#     Cfc18 = zeros(n)
#     Cfc19 = zeros(n)

#     for i=1:n
#         ti = tvec[i]

#         ### Calculate constants
#         M = V(ti)/a #Mach number
#         beta = sqrt(1-(M^2)) #Compressibility factor

#         ### Calculate time constants
#         Tl = c/a
#         k_alpha = 1/((1-M) + (pi*beta*(M^2)*(A[1]*b[1] + A[2]*b[2]))) #I wonder if any of these should include slope. -> Nope, they are included in the normal force equations. 

#         k_q = 1/((1-M) + (2*pi*beta*(M^2)*(A[1]*b[1] + A[2]*b[2])))

#         ### Attached flow
#         #Circulatory normal force
#         Ccn[i] = dcldalpha*(2*V(ti)/c)*(beta^2)*(A[1]*b[1]*u[i,1] + A[2]*b[2]*u[i,2])

#         #noncirculatory normal force due to angle of attack
#         Cina[i] = 4*du[i,3]/M 

#         #noncirculatory normal force due to pitch rate
#         Cinq[i] = du[i,4]/M 

#         #Total attached flow normal force
#         Cpn[i] = Ccn[i] + Cina[i] + Cinq[i]

#         ### Dynamic Stall  
#         alphaE[i] = (beta^2)*(2*V(ti)/c)*(A[1]*b[1]*u[i,1] + A[2]*b[2]*u[i,2]) #Equivalent angle of attack (due to shed wake terms)
#         K[i] = ((1+sqrt(u[i,6]))^2)/4
#         Cf[i] = dcldalpha*K[i]*alphaE[i] + Cina[i] + Cinq[i] 
#         Cdn[i] = Cf[i] + u[i,8] #Total normal force - dynamic stall conditions

#         ### Drag
#         alpha_f = u[i,5]/dcldalpha
#         ff = separationpoint_BL(alpha_f, alpha1, alpha0)
#         f = separationpoint_BL(alpha(ti), alpha1, alpha0)
#         Cc[i] = eta*dcldalpha*(alphaE[i]^2)*sqrt(f)
#         Cfc[i]  = eta*dcldalpha*((alphaE[i])^2)*sqrt(u[i,6])
#         Cfc2[i] = eta*dcldalpha*((alphaE[i])^2)*sqrt(ff) 
#         Cfc3[i] = eta*alphaE[i]*(dcldalpha*K[i]*alphaE[i] + Cina[i] + Cinq[i])*sqrt(u[i,6]) #DSF
#         Cfc4[i] = eta*tan(alphaE[i])*(Cf[i])
#         Cfc5[i] = eta*Cdn[i]*sin(alphaE[i])*cos(alphaE[i])
#         Cfc6[i] = Ccn[i]*eta*sqrt(u[i,6])*cos(alphaE[i])
#         Cfc7[i] = eta*dcldalpha*sin(alphaE[i])*alphaE[i]*(sqrt(f))
#         Cfc8[i] = eta*dcldalpha*sin(alphaE[i])*alphaE[i]*sqrt(u[i,6])
#         Cfc9[i] = eta*sin(alphaE[i])*Cf[i]*sqrt(u[i,6])*cos(alpha(ti)) #DSF
#         Cfc10[i] = eta*dcldalpha*((alphaE[i])^2)*(1-sqrt(u[i,6])) 
#         Cfc11[i] = eta*alphaE[i]*Cdn[i]*sqrt(u[i,6])
#         Cfc12[i] = eta*Cdn[i]*sin(alphaE[i])/cos(alpha(ti))
#         Cfc13[i] = eta*Cdn[i]*sin(alphaE[i])/cos(alphaE[i])
#         Cfc14[i] = eta*Ccn[i]*sin(alphaE[i])/cos(alpha(ti))
#         Cfc15[i] = eta*sqrt(u[i,6])*Cdn[i]*sin(alphaE[i])/cos(alpha(ti))
#         Cfc16[i] = eta*K[i]*Cdn[i]*sin(alphaE[i])/cos(alpha(ti)) 
#         Cfc17[i] = eta*Cf[i]*sin(alphaE[i])/cos(alpha(ti))
#         Cfc18[i] = eta*sqrt(u[i,6])*Cf[i]*sin(alphaE[i])/cos(alpha(ti)) #double includes separation factor. 
#         Cfc19[i] = eta*sqrt(u[i,6])*Ccn[i]*sin(alphaE[i])/cos(alpha(ti))

#         Cdc[i] = Cfc[i] #It think there should be a vortex component.
#         # println("")
#         # println("time step: ", i)
#         # println(" cos(alpha): ", cos(alpha(ti)))
#         # println("cos(alphaE): ", cos(alphaE[i]))
#         # println("alpha(t): ", alpha(ti))
#         # println("  alphaE: ", alphaE[i]) #AlphaE lags behind alpha

#         Cdd[i] = Cdn[i]*sin(alpha(ti)) - Cdc[i]*cos(alpha(ti))
#         Cdd2[i] = Cdc[i]*cos(alpha(ti))
#         # Cdd3[i] = -Cdc[i]*cos(alphaE[i])*log(alpha(ti)) #(alpha(ti)^(1/3))
#         Cdd4[i] = Cfc2[i]*cos(alpha(ti))
#         Cdd5[i] = Cdn[i]*sin(alpha(ti)) - Cfc2[i]*cos(alpha(ti))
#         Cdd6[i] = Cdn[i]*sin(alpha(ti))
#         Cdd7[i] = Cdn[i]*sin(alpha(ti)) - Cfc3[i]*cos(alpha(ti))
#         Cdd8[i] = Cdn[i]*sin(alphaE[i]) - Cfc3[i]*cos(alphaE[i])
#         Cdd9[i] = Cdn[i]*sin(alpha(ti)) - Cfc4[i]*cos(alpha(ti))
#         Cdd10[i] = Cfc5[i]*cos(alpha(ti))
#         Cdd11[i] = Cdn[i]*sin(alpha(ti)) - Cfc7[i]*cos(alpha(ti))
#         Cdd12[i] = Cdn[i]*sin(alpha(ti)) - Cfc8[i]*cos(alpha(ti))
#         Cdd13[i] = Cfc10[i]*cos(alphaE[i])
#         Cdd14[i] = Cdn[i]*sin(alpha(ti)) - Cfc11[i]*cos(alpha(ti))
#         Cdd15[i] = Cdn[i]*sin(alphaE[i])
#         Cdd16[i] = Cdn[i]*sin(alpha(ti)) - Cfc12[i]*cos(alpha(ti))
#         Cdd17[i] = Cdn[i]*sin(alpha(ti)) - Cfc13[i]*cos(alpha(ti))
#         Cdd18[i] = Cdn[i]*sin(alpha(ti)) - Cfc14[i]*cos(alpha(ti))
#         Cdd19[i] = Cdn[i]*sin(alpha(ti)) - Cfc19[i]*cos(alpha(ti))
#     end
#     return Cdn, Cpn, Cdd, Cdd2, Cdd3, Cdd4, Cdd5, Cdd6, Cdd7, Cdd8, Cdd9, Cdd10, Cdd11, Cdd12, Cdd13, Cdd14, Cdd15, Cdd16, Cdd17, Cdd18, Cdd19, u, du, Cc, Cfc, Cfc2, Cfc3, Cfc4, Cfc5, Cfc6, Cfc7, Cfc8, Cfc9, Cfc10, Cfc11, Cfc12, Cfc13, Cfc14, Cfc15, Cfc16, Cfc17, Cfc18, Cfc19
# end

function parsesolution_states(sol, p; eta=0.95)
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

function fullstates_parsesolution(sol, p; eta=0.95)
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

    V, Vdot, alpha, alphadot, c, dCndalpha, alpha1, Cn1, Cmo, A, b, S, K, m, T_p, T_f, T_vl, T_v, a = p

    t, x, dx = extractdata(sol)

    ### Constants
    L = 2*V(t)/c
    tau_p = T_p/L
    tau_f = T_f/L
    tau_vl = T_vl/L
    tau_v = T_v/L

    n = length(t)

    #Initialize arrays
    Ccn = zeros(n)
    Cina = zeros(n)
    Cinq = zeros(n)
    Cpn = zeros(n)

    alphaE = zeros(n)
    Cfn = zeros(n)
    Cfc = zeros(n)
    Cfm = zeros(n)
    Cvm = zeros(n)

    Cnt = zeros(n)
    Cnt2 = zeros(n)
    Cdt = zeros(n)
    Cmt = zeros(n)

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

        Cfn[i] = dCndalpha*(((1+sqrt(x[i,10]))^2)/4)*alphaE[i] + Cina[i] + Cinq[i] #Nonlinear dynamic normal force
        Cfc[i] = eta*dCndalpha*(alphaE[i]^2)*sqrt(x[i,10]) #Chordwise Force
        Cfm[i] = (K[1] + (K[2]*(1-x[i,13])) + (K[3]*sin(pi*(x[i,13]^m))))*Ccn[i] + Cmo

        # Cpv = (1 - cos(pi*x[i,11]/tau_vl))/4
        Cpv = (1 - cos(pi*x[i,11]/T_vl))/4
        Cvm[i] = -Cpv*x[i,12]

        Cnt[i] = Cfn[i] + x[i,12] #Total dynamic normal force
        Cnt2[i] = (((1+sqrt(x[i,10]))^2)/4)*(dCndalpha*alphaE[i] + Cina[i] + Cinq[i]) + x[i,12] #As above but including the nonlinear forces in the diminished portion of the lift. 
        Cdt[i] = Cnt[i]*sin(alpha(t[i])) - Cfc[i]*cos(alpha(t[i])) #Drag
        Cmt[i] = Cfm[i] + Cvm[i]
    end
    return t, x, dx, Cnt, Cdt, Cmt
end



function analyzesolution(sol, p; eta=0.95)
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

    V, Vdot, alpha, alphadot, c, dCndalpha, alpha1, Cn1, Cmo, A, b, S, K, m, T_p, T_f, T_vl, T_v, a = p

    t, x, dx = extractdata(sol)

    ### Constants
    L = 2*V(t)/c
    tau_p = T_p/L
    tau_f = T_f/L
    tau_vl = T_vl/L
    tau_v = T_v/L

    n = length(t)

    #Initialize arrays
    Ccn = zeros(n)
    Cina = zeros(n)
    Cinq = zeros(n)
    Cpn = zeros(n)

    Ccndot = zeros(n)
    Cvdot = zeros(n)
    Cv = zeros(n)
    Hcv = zeros(n)

    alphaE = zeros(n)
    Cfn = zeros(n)
    Cfc = zeros(n)
    Cfm = zeros(n)
    Cvm = zeros(n)

    Cnt = zeros(n)
    Cnt2 = zeros(n)
    Cdt = zeros(n)
    Cmt = zeros(n)

    tauflag = true
    tausflag = false
    tvlflag = true
    tstartvec = []
    tstopvec = []
    tvlvec = []

    for i = 1:n
        ### Calculate constants
        M = V(t[i])/a
        beta = sqrt(1-(M^2))

        ### Calculate attached lift
        Ccn[i] = dCndalpha*(2*V(t[i])/c)*(beta^2)*((A[1]*b[1]*x[i,1])+(A[2]*b[2]*x[i,2])) #Circulatory lift
        Cina[i] = (4/M)*dx[i,3] #Normal force due to angle of attack  
        Cinq[i] = (1/M)*dx[i,4] #Normal force due to pitch rate 

        Cpn[i] = Ccn[i] + Cina[i] + Cinq[i] #Attached normal force

        Ccndot[i] = dCndalpha*(2*Vdot(t)/c)*(1-(M^2))*((A[1]*b[1]*x[i,1])+(A[2]*b[2]*x[i,2])) + dCndalpha*(2*V(t)/c)*(1-(M^2))*((A[1]*b[1]*dx[i,1])+(A[2]*b[2]*dx[i,2])) - dCndalpha*(4*(M^2)/c)*Vdot(t)*((A[1]*b[1]*x[i,1])+(A[2]*b[2]*x[i,2])) #Derivative of the circulatory normal force with respect to time. #This matches the numerical derivative of Ccn. 

        # Cvdot = (Ccndot*(1-((1+sqrt(x[10]))^2)/4) - Ccn*(1 + sqrt(x[10]))*dx[10]/(4*sqrt(x[10])))*Heavi(2*tau_vl - x[11])*Heavi(x[9]-Cn1) #*(c/a) #*(2*V(t)/c) #Derivative of the vortex lift contribution with respect to time. #Todo: There is something going on here. I don't think the units quite work out. I tried correcting the units by multiplying by 2V/c, but that gave me this huge jump. And I don't think there is supposed to be this huge jump. At least, that's what the paper's solutions show. Thus there is something up. I also tried multiplying by c/a... which doesn't help because it makes the forcing function even smaller. 
        # Cvdot = Ccn*(1-(((1+sqrt(x[10]))/2)^2))*Heavi(2*tau_vl - x[11])*Heavi(x[9]-Cn1) #Note: Cvn dropped even more. 
        # Cvdot[i] = (Ccndot[i]*(1-((1+sqrt(x[i,10]))^2)/4) - Ccn[i]*(1 + sqrt(x[i,10]))*dx[i,10]/(4*sqrt(x[i,10])))*Heavi(2*tau_vl - x[i,11]) - Ccn[i]*(1-((1+sqrt(x[i,10]))^2)/4)*dx[i,11]*diracdelta(2*tau_vl-x[i,11])
        # Cvdot[i] = (Ccndot[i]*(1-((1+sqrt(x[i,10]))^2)/4) - Ccn[i]*(1 + sqrt(x[i,10]))*dx[i,10]/(4*sqrt(x[i,10])))*Heavi(2*tau_vl - x[i,11])
        Cvdot[i] = (Ccndot[i]*(1-((1+sqrt(x[i,10]))^2)/4) - Ccn[i]*(1 + sqrt(x[i,10]))*dx[i,10]/(4*sqrt(x[i,10])))*Heavi(2*T_vl - x[i,11]) #Note: This was the one. 
        # Cvdot[i] = (Ccndot[i]*(1-((1+sqrt(x[i,10]))^2)/4) - Ccn[i]*(1 + sqrt(x[i,10]))*dx[i,10]/(4*sqrt(x[i,10])))*Heavi(T_vl - x[i,11]) - Ccn[i]*(1-((1+sqrt(x[i,10]))^2)/4)*dx[i,11]*dirac(T_vl-x[i,11])

        # Cv[i] = Ccn[i]*(1-((1+sqrt(x[i,10]))^2)/4)*Heavi(2*tau_vl - x[i,11])
        Cv[i] = Ccn[i]*(1-((1+sqrt(x[i,10]))^2)/4)*Heavi(2*T_vl - x[i,11])

        # Hcv[i] = Heavi(2*tau_vl - x[i,11])
        Hcv[i] = Heavi(2*T_vl - x[i,11])


        alphaE[i] = (beta^2)*(2*V(t[i])/c)*((A[1]*b[1]*x[i,1])+(A[2]*b[2]*x[i,2])) #Effective angle of attack

        Cfn[i] = dCndalpha*(((1+sqrt(x[i,10]))^2)/4)*alphaE[i] + Cina[i] + Cinq[i] #Nonlinear dynamic normal force
        Cfc[i] = eta*dCndalpha*(alphaE[i]^2)*sqrt(x[i,10]) #Chordwise Force
        Cfm[i] = (K[1] + (K[2]*(1-x[i,13])) + (K[3]*sin(pi*(x[i,13]^m))))*Ccn[i] + Cmo

        # Cpv = (1 - cos(pi*x[i,11]/tau_vl))/4
        Cpv = (1 - cos(pi*x[i,11]/T_vl))/4
        Cvm[i] = -Cpv*x[i,12]

        Cnt[i] = Cfn[i] + x[i,12] #Total dynamic normal force
        Cnt2[i] = (((1+sqrt(x[i,10]))^2)/4)*(dCndalpha*alphaE[i] + Cina[i] + Cinq[i]) + x[i,12] #As above but including the nonlinear forces in the diminished portion of the lift. 
        Cdt[i] = Cnt[i]*sin(alpha(t[i])) - Cfc[i]*cos(alpha(t[i])) #Drag
        Cmt[i] = Cfm[i] + Cvm[i]

        if x[i,11]>=2*T_vl && tvlflag
            tvlflag = false
            push!(tvlvec, t[i])
        end
        if x[i,11]>0 && tauflag
            tauflag = false
            tausflag = true
            push!(tstartvec, t[i])
        elseif isapprox(x[i,11],0.0, atol=0.0001) && tausflag
            tauflag = true
            tausflag = false
            tvlflag = true
            # println("Got here")
            push!(tstopvec, t[i])
        end
    end
    return Cv, Cvdot, Ccn, Ccndot, tstartvec, tstopvec, tvlvec, Hcv
end


