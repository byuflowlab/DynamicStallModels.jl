using FLOWMath

function Heavi(t)
    if t>0
        return 1
    else
        return 0
    end
end

function ffun(alphan, alpha1, alpha0) #We could try some other f functions.
    if alphan<=alpha1
        return 1.0 - 0.3*exp((alphan-alpha1)/S[1])
    else
        return 0.04 + 0.66*exp((alpha1-alphan)/S[2])
    end
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
    a33 = -1/(k_alpha*Tl)
    a44 = -1/(k_q*Tl)
    a55 = -1/(b[3]*k_alpha_m*Tl)
    a66 = -1/(b[4]*k_alpha_m*Tl)
    a77 = -2*b[5]*V*(beta^2)/c
    a88 = -1/(k_q_m*Tl)

    ### Calculate the state rates 
    # I decided to leave the state rates as a 1D array rather than a matrix. A matrix would have been easier to code... but I'm not sure if it would have been more computationally efficient. 
    du[1] = u[1]*a11 + alpha(t) + q(t)/2
    du[2] = u[2]*a22 + alpha(t) + q(t)/2
    du[3] = u[3]*a33 + alpha(t)
    du[4] = u[4]*a44 + q(t)
    du[5] = u[5]*a55 + alpha(t)
    du[6] = u[6]*a66 + alpha(t)
    du[7] = u[7]*a77 + q(t)
    du[8] = u[8]*a88 + q(t)

end

function states_liftonly!(du, u, p, t)
    ### Unpack p
    V, alpha, alphadot, c, dcldalpha, alpha1, alpha0, Cn1, A, b, a, Tp, Tf, Tv, Tvl = p

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
    a33 = -1/(k_alpha*Tl)
    a44 = -1/(k_q*Tl)
    
    ### Calculate the attached state rates 
    #Shed wake terms
    du[1] = u[1]*a11 + alpha(t) + q(t)/2
    du[2] = u[2]*a22 + alpha(t) + q(t)/2

    #noncirculatory angle of attack contribution
    du[3] = u[3]*a33 + alpha(t)

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
    du[5] = -u[5]/Tp + Cpn/Tp #sv9 #(Cpn-u[5])/Tp
    
    # Delayed effective angle of attack
    alpha_f = u[5]/dcldalpha

    # Unsteady TE separation point, f''
    du[6] = -u[6]/Tf + ffun(alpha_f, alpha1, alpha0)/Tf #sv10 

    tau_vl = 2*Tvl

    # Separation point location time constant, τ
    du[7] = V(t)*Heavi(u[5]-Cn1)*Heavi(tau_vl-u[7])/(3*c) - V(t)*Heavi(Cn1-u[5])*Heavi(u[7])/c #Note: This is my own creation. Especially the stuff after the negative sign. I figure it will take some time to reattach fully, but I do know that there can be multiple vortices attached... so it probably should just reset to 0. Or maybe reset faster. #TODO: I need to check that the Heavi function is performing as I want it to. 

    # Change in lift due to vortex position
    Ccndot = dcldalpha*(2*V(t)/c)*(beta^2)*(A[1]*b[1]*du[1] + A[2]*b[2]*du[2]) 
    Cvdot = (Ccndot*(1- ((1+sqrt(u[6]))^2)/4) - Ccn*(1+sqrt(u[6]))*du[6]/(4*sqrt(u[6])))*Heavi(tau_vl-u[7])*Heavi(u[5]-Cn1) #Note: If u[6]<=0, then we be screwed. (u[6]=f'') 

    # Vortex lift contribution
    du[8] = -u[8]/Tv + Cvdot/Tv  
end

function parsesolution_liftonly(sol, p)
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

    Cl = zeros(n)
    Cd = zeros(n)

    Cn = zeros(n)
    Cc = zeros(n)

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

    for i=1:n
        ti = tvec[i]

        ### Calculate constants
        M = V(ti)/a #Mach number
        beta = sqrt(1-(M^2)) #Compressibility factor

        ### Calculate time constants
        Tl = c/a
        k_alpha = 1/((1-M) + (pi*beta*(M^2)*(A[1]*b[1] + A[2]*b[2]))) #I wonder if any of these should include slope. -> Nope, they are included in the normal force equations. 

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

    end
    return Cdn
end

function parsesolution_liftonly(sol, p)
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

    Cl = zeros(n)
    Cd = zeros(n)

    Cn = zeros(n)
    Cc = zeros(n)

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

    for i=1:n
        ti = tvec[i]

        ### Calculate constants
        M = V(ti)/a #Mach number
        beta = sqrt(1-(M^2)) #Compressibility factor

        ### Calculate time constants
        Tl = c/a
        k_alpha = 1/((1-M) + (pi*beta*(M^2)*(A[1]*b[1] + A[2]*b[2]))) #I wonder if any of these should include slope. -> Nope, they are included in the normal force equations. 

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
    end
    return Cdn, Cpn, u, du
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

    q(t) = alphadot(t)*c/V(t)

    V, alpha, alphadot, c, dcldalpha, alpha1, alpha0, Cn1, A, b, a, Tp, Tf, Tv = p

    tvec, u, du = extractdata(sol)
    n = length(tvec)

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

    Cc = zeros(n)
    Cdc = zeros(n)

    Cdd = zeros(n)
    Cdd2 = zeros(n)
    Cdd3 = zeros(n)
    Cdd4 = zeros(n)
    Cdd5 = zeros(n)
    Cdd6 = zeros(n)
    Cdd7 = zeros(n)
    Cdd8 = zeros(n)
    Cdd9 = zeros(n)
    Cdd10 = zeros(n)
    Cdd11 = zeros(n)
    Cdd12 = zeros(n)
    Cdd13 = zeros(n)
    Cdd14 = zeros(n)
    Cdd15 = zeros(n)
    Cdd16 = zeros(n)
    Cdd17 = zeros(n)
    Cdd18 = zeros(n)
    Cdd19 = zeros(n)

    Cfc = zeros(n)
    Cfc2 = zeros(n)
    Cfc3 = zeros(n)
    Cfc4 = zeros(n)
    Cfc5 = zeros(n)
    Cfc6 = zeros(n)
    Cfc7 = zeros(n)
    Cfc8 = zeros(n)
    Cfc9 = zeros(n)
    Cfc10 = zeros(n)
    Cfc11 = zeros(n)
    Cfc12 = zeros(n)
    Cfc13 = zeros(n)
    Cfc14 = zeros(n)
    Cfc15 = zeros(n)
    Cfc16 = zeros(n)
    Cfc17 = zeros(n)
    Cfc18 = zeros(n)
    Cfc19 = zeros(n)

    for i=1:n
        ti = tvec[i]

        ### Calculate constants
        M = V(ti)/a #Mach number
        beta = sqrt(1-(M^2)) #Compressibility factor

        ### Calculate time constants
        Tl = c/a
        k_alpha = 1/((1-M) + (pi*beta*(M^2)*(A[1]*b[1] + A[2]*b[2]))) #I wonder if any of these should include slope. -> Nope, they are included in the normal force equations. 

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
        alpha_f = u[i,5]/dcldalpha
        ff = ffun(alpha_f, alpha1, alpha0)
        f = ffun(alpha(ti), alpha1, alpha0)
        Cc[i] = eta*dcldalpha*(alphaE[i]^2)*sqrt(f)
        Cfc[i]  = eta*dcldalpha*((alphaE[i])^2)*sqrt(u[i,6])
        Cfc2[i] = eta*dcldalpha*((alphaE[i])^2)*sqrt(ff) 
        Cfc3[i] = eta*alphaE[i]*(dcldalpha*K[i]*alphaE[i] + Cina[i] + Cinq[i])*sqrt(u[i,6]) #DSF
        Cfc4[i] = eta*tan(alphaE[i])*(Cf[i])
        Cfc5[i] = eta*Cdn[i]*sin(alphaE[i])*cos(alphaE[i])
        Cfc6[i] = Ccn[i]*eta*sqrt(u[i,6])*cos(alphaE[i])
        Cfc7[i] = eta*dcldalpha*sin(alphaE[i])*alphaE[i]*(sqrt(f))
        Cfc8[i] = eta*dcldalpha*sin(alphaE[i])*alphaE[i]*sqrt(u[i,6])
        Cfc9[i] = eta*sin(alphaE[i])*Cf[i]*sqrt(u[i,6])*cos(alpha(ti)) #DSF
        Cfc10[i] = eta*dcldalpha*((alphaE[i])^2)*(1-sqrt(u[i,6])) 
        Cfc11[i] = eta*alphaE[i]*Cdn[i]*sqrt(u[i,6])
        Cfc12[i] = eta*Cdn[i]*sin(alphaE[i])/cos(alpha(ti))
        Cfc13[i] = eta*Cdn[i]*sin(alphaE[i])/cos(alphaE[i])
        Cfc14[i] = eta*Ccn[i]*sin(alphaE[i])/cos(alpha(ti))
        Cfc15[i] = eta*sqrt(u[i,6])*Cdn[i]*sin(alphaE[i])/cos(alpha(ti))
        Cfc16[i] = eta*K[i]*Cdn[i]*sin(alphaE[i])/cos(alpha(ti)) 
        Cfc17[i] = eta*Cf[i]*sin(alphaE[i])/cos(alpha(ti))
        Cfc18[i] = eta*sqrt(u[i,6])*Cf[i]*sin(alphaE[i])/cos(alpha(ti)) #double includes separation factor. 
        Cfc19[i] = eta*sqrt(u[i,6])*Ccn[i]*sin(alphaE[i])/cos(alpha(ti))

        Cdc[i] = Cfc[i] #It think there should be a vortex component.
        # println("")
        # println("time step: ", i)
        # println(" cos(alpha): ", cos(alpha(ti)))
        # println("cos(alphaE): ", cos(alphaE[i]))
        # println("alpha(t): ", alpha(ti))
        # println("  alphaE: ", alphaE[i]) #AlphaE lags behind alpha

        Cdd[i] = Cdn[i]*sin(alpha(ti)) - Cdc[i]*cos(alpha(ti))
        Cdd2[i] = Cdc[i]*cos(alpha(ti))
        # Cdd3[i] = -Cdc[i]*cos(alphaE[i])*log(alpha(ti)) #(alpha(ti)^(1/3))
        Cdd4[i] = Cfc2[i]*cos(alpha(ti))
        Cdd5[i] = Cdn[i]*sin(alpha(ti)) - Cfc2[i]*cos(alpha(ti))
        Cdd6[i] = Cdn[i]*sin(alpha(ti))
        Cdd7[i] = Cdn[i]*sin(alpha(ti)) - Cfc3[i]*cos(alpha(ti))
        Cdd8[i] = Cdn[i]*sin(alphaE[i]) - Cfc3[i]*cos(alphaE[i])
        Cdd9[i] = Cdn[i]*sin(alpha(ti)) - Cfc4[i]*cos(alpha(ti))
        Cdd10[i] = Cfc5[i]*cos(alpha(ti))
        Cdd11[i] = Cdn[i]*sin(alpha(ti)) - Cfc7[i]*cos(alpha(ti))
        Cdd12[i] = Cdn[i]*sin(alpha(ti)) - Cfc8[i]*cos(alpha(ti))
        Cdd13[i] = Cfc10[i]*cos(alphaE[i])
        Cdd14[i] = Cdn[i]*sin(alpha(ti)) - Cfc11[i]*cos(alpha(ti))
        Cdd15[i] = Cdn[i]*sin(alphaE[i])
        Cdd16[i] = Cdn[i]*sin(alpha(ti)) - Cfc12[i]*cos(alpha(ti))
        Cdd17[i] = Cdn[i]*sin(alpha(ti)) - Cfc13[i]*cos(alpha(ti))
        Cdd18[i] = Cdn[i]*sin(alpha(ti)) - Cfc14[i]*cos(alpha(ti))
        Cdd19[i] = Cdn[i]*sin(alpha(ti)) - Cfc19[i]*cos(alpha(ti))
    end
    return Cdn, Cpn, Cdd, Cdd2, Cdd3, Cdd4, Cdd5, Cdd6, Cdd7, Cdd8, Cdd9, Cdd10, Cdd11, Cdd12, Cdd13, Cdd14, Cdd15, Cdd16, Cdd17, Cdd18, Cdd19, u, du, Cc, Cfc, Cfc2, Cfc3, Cfc4, Cfc5, Cfc6, Cfc7, Cfc8, Cfc9, Cfc10, Cfc11, Cfc12, Cfc13, Cfc14, Cfc15, Cfc16, Cfc17, Cfc18, Cfc19
end