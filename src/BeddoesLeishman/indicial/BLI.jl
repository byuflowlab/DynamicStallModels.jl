function Heavi(t)
    if t>0
        return 1
    else
        return 0
    end
end

function ffun(alphan, alpha1, S) #x
    if alphan<=alpha1
        return 1.0 - 0.3*exp((alphan-alpha1)/S[1])
    else
        return 0.04 + 0.66*exp((alpha1-alphan)/S[2])
    end
end

function BeddoesLeishman(alphafun, ufun, tspan, c, A, b, S, dCndalpha, Cn1, alpha1, Tf, Tp, Tvl; n=100, a=343.3, f0=NaN)

    t = collect(range(tspan[1], tspan[2], length=n))
    alpha = alphafun.(t) #Todo: This isn't oscillating.
    v = ufun.(t)
    Tl = c/a

    delS = zeros(n)
    delalpha = zeros(n)
    delt = zeros(n)
    M = zeros(n)
    beta = zeros(n)
    X = zeros(n)
    Y = zeros(n)
    alphaE = zeros(n)
    alphaf = zeros(n)
    Ccn = zeros(n)
    kalpha = zeros(n)
    D = zeros(n)
    Cin = zeros(n)
    Cpn = zeros(n)
    fp = zeros(n)
    Df = zeros(n)
    fpp = zeros(n)
    Cfn = zeros(n)
    Dp = zeros(n)
    Cprimen = zeros(n)
    k = zeros(n)
    tauv = zeros(n)
    Cv = zeros(n)
    Cvn = zeros(n)
    Cn = zeros(n)

    ### Run initial time step
    delS[1] = 2*(v[1]*t[1] - v[2]*t[2])/c #Distance traveled in semi-chords
    delalpha[1] = alpha[1] - alpha[2]
    delt[1] = t[1] - t[2]
    M[1] = v[1]/a #Mach number
    beta[1] = sqrt(1 - M[1]^2) #Compressibility factor

    X[1] = A[1]*delalpha[1]*exp(-b[1]*(beta[1]^2)*delS[1]/2)
    Y[1] = A[2]*delalpha[1]*exp(-b[2]*(beta[1]^2)*delS[1]/2)

    alphaE[1] = alpha[1] - X[1] - Y[1] #Equivalent Angle of Attack
    Ccn[1] = dCndalpha*alphaE[1]

    kalpha[1] = 0.75/((1-M[1]) + (pi*(beta[1]^2)*(M[1]^2)*(A[1]*b[1] + A[2]*b[2])))
    D[1] = ((delalpha[1] - delalpha[2])/delt[1])*exp(-delt[1]/(2*kalpha[1]*Tl)) #nonlinear defeciency function

    Cin[1] = 4*kalpha[1]*Tl*((delalpha[1]/delt[1]) - D[1])/M[1] #Nonlinear normal force contribution
    Cpn[1] = Ccn[1] + Cin[1]

    Cprimen[1] = Cpn[1]

    alphaf[1] = alpha[1]
    if isnan(f0)
        fp[1] = ffun(alphaf[1], alpha1, S) #Attachment factor
    else
        fp[1] = f0
    end
    Df[1] = (1 - fp[1])*exp(delS[1]/(2*Tf)) #Attachment factor defeciency function
    fpp[1] = fp[1] - Df[1] #Delayed attachment factor
    Cfn[1] = dCndalpha*((1 + sqrt(fpp[1]))^2)*alphaE[1]/4 + Cin[1] #Nonlinear dynamic normal force

    k[1] = ((1 + sqrt(fpp[1]))^2)/4
    tauv[1] = 0.0
    Cv[1] = Ccn[1]*(1-k[1])

    Cvn[1] = 0.0

    Cn[1] = Cfn[1] + Cvn[1]

    ### Run the rest of the time steps
    for i=2:n
        delS[i] = 2*(v[i]*t[i] - v[i-1]*t[i-1])/c #Distance traveled in semi-chords
        delalpha[i] = alpha[i] - alpha[i-1]
        delt[i] = t[i] - t[i-1]
        M[i] = v[i]/a #Mach number
        beta[i] = sqrt(1.0-M[i]^2) #Compressibility factor

        X[i] = X[i-1]*exp(-b[1]*(beta[i]^2)*delS[i]) + A[1]*delalpha[i]*exp(-b[1]*(beta[i]^2)*delS[i]/2)
        Y[i] = Y[i-1]*exp(-b[2]*(beta[i]^2)*delS[i]) + A[2]*delalpha[i]*exp(-b[2]*(beta[i]^2)*delS[i]/2)

        alphaE[i] = alpha[i] - X[i] - Y[i] #Equivalent Angle of Attack
        # println("")
        # println("t-maximin: ", extrema(t))
        # println("maximin: ", extrema(alpha))
        # println("αn: ", alpha[i])
        # println("X: ", X[i])
        # println("Y: ", Y[i])

        Ccn[i] = dCndalpha*alphaE[i] #Noraml force due to Circulatory effects  #Todo: This is basically the same throughout. 

        kalpha[i] = 0.75/((1-M[i]) + (pi*(beta[i]^2)*(M[i]^2)*(A[1]*b[1] + A[2]*b[2])))
        D[i] = D[i-1]*exp(-delt[i]/(kalpha[i]*Tl)) + ((delalpha[i] - delalpha[i-1])/delt[i])*exp(-delt[i]/(2*kalpha[i]*Tl)) #nonlinear defeciency function
        Cin[i] = 4*kalpha[i]*Tl*((delalpha[i]/delt[i]) - D[i])/M[i] #Nonlinear normal force contribution

        Cpn[i] = Ccn[i] + Cin[i]

        Dp[i] = Dp[i-1]*exp(-delS[i]/Tp) + (Cpn[i]-Cpn[i-1])*exp(-delS[i]/(2*Tp))
        Cprimen[i] = Cpn[i] - Dp[i]

        alphaf[i] = Cprimen[i]/dCndalpha
        fp[i] = ffun(alphaf[i], alpha1, S) #Attachment factor Todo: alphaf???
        Df[i] = Df[i-1]*exp(-delS[i]/Tf) + (fp[i] - fp[i-1])*exp(-delS[i]/(2*Tf)) #Attachment factor defeciency function #Note: The paper doesn't have a negative in the exponent, but other papers do. Which makes sense, for a decaying first order ODE. 
        fpp[i] = fp[i] - Df[i] #Delayed attachment factor

        Cfn[i] = dCndalpha*((1 + sqrt(fpp[i]))^2)*alphaE[i]/4 + Cin[i] #Nonlinear dynamic normal force #Todo: This has essentially the same value throughout. 

        k[i] = ((1 + sqrt(fpp[i]))^2)/4

        # tauv[i] = (tauv[i-1] + delt[i]*v[i]*Heavi(Tvl-tauv[i])/(2*c))*Heavi(Cprimen[i]-Cn1) #Todo. I'm not sure if this is right. I don't think there should be two heaviside functions in here. 
        tauv[i] = (tauv[i-1] + delt[i]*v[i])*Heavi(Cprimen[i]-Cn1)
        Cv[i] = Ccn[i]*(1-k[i])
        Cvn[i] = Cvn[i-1]*exp(-delS[i]/Tvl) + (Cv[i] - Cv[i-1])*exp(-delS[i]/(2*Tvl))*Heavi(Tvl-tauv[i]) #Todo: This is super small. But that's expected, since this appears to be attached conditions. 

        Cn[i] = Cfn[i] + Cvn[i]
    end
    # println(Cvn)
    return Cn, t, Cpn, Ccn, Cin
end