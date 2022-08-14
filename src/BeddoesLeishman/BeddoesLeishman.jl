using FLOWMath

#=
5/5/2022 

I found this in a "semester" directory. It needs to be gone through and compared to my current work. 



=#







function Heavi(t)
    if t>0
        return 1
    else
        return 0
    end
end

function ffun(alphan, alpha1, alpha0)
    if alphan<=alpha1
        return 1.0 - 0.3*exp((alphan-alpha1)/S[1])
    else
        return 0.04 + 0.66*exp((alpha1-alphan)/S[2])
    end
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
    du[5] = (Cpn-u[5])/Tp #-u[5]/Tp + Cpn/Tp #sv9
    
    # Delayed effective angle of attack
    alpha_f = u[5]/dcldalpha

    # Unsteady TE separation point, f''
    du[6] = -u[6]/Tf + ffun(alpha_f, alpha1, alpha0)/Tf #sv10 

    tau_v = 2*Tv

    # Separation point location time constant, τ
    du[7] = V(t)*Heavi(u[5]-Cn1)*Heavi(tau_v-u[7])/(3*c) - V(t)*Heavi(Cn1-u[5])*Heavi(u[7])/c #Note: This is my own creation. Especially the stuff after the negative sign. I figure it will take some time to reattach fully, but I do know that there can be multiple vortices attached... so it probably should just reset to 0. Or maybe reset faster.  

    # Change in lift due to vortex position
    Ccndot = dcldalpha*(2*V(t)/c)*(beta^2)*(A[1]*b[1]*du[1] + A[2]*b[2]*du[2]) 
    Cvdot = (Ccndot*(1- ((1+sqrt(u[6]))^2)/4) - Ccn*(1+sqrt(u[6]))*du[6]/(4*sqrt(u[6])))*Heavi(tau_v-u[7])*Heavi(u[5]-Cn1) #Note: If u[6]<=0, then we be screwed. (u[6]=f'') 

    # Vortex lift contribution
    du[8] = -u[8]/Tv + Cvdot/Tv  
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