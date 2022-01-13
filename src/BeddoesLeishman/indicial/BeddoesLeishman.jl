



# function alpha1fun() #Not really sure what to do here... cause f is a function of alpha and alpha1... and I'm pretty sure alpha1 is constant for an airfoil... if I remember correctly. => They called it the static stall angle, so it's probably just user input of when we it stalls. 
#     ffun(alpha)
# end

function ffun(alphan, alpha1, S1, S2) #x
    if alphan<=alpha1
        return 1.0 - 0.3*exp((alphan-alpha1)/S1)
    else
        return 0.04 + 0.66*exp((alpha1-alphan)/S2)
    end
end

function KNfun(fppn) #x
    return ((1 + sqrt(fppn))^2)/4
end

function Cvfun(KNn, CcNn) #x
    return CcNn*(1-KNn)
end

function CvNfun(CvNn_1, Cvn, Cvn_1, delS, Tv) #x
    return CvNn_1*exp(delS/Tv) + (Cvn-Cvn_1)*exp(delS/(2*Tv))
end

function Dffun(Dfn_1, delSn, Tf, fpn, fpn_1) #x
    return Dfn_1*exp(delSn/Tf) + (fpn-fpn_1)*exp(delSn/(2*Tf))
end

function delSfun(c, Un, Un_1, tn, tn_1) #x #Note: I recalculated this myself, but I thought I saw a different formulation of this somewhere. 
    return 2*((Un*tn)-(Un_1*tn_1))/c
end

function betafun(Mn) #x
    return sqrt(1-(Mn^2))
end

function Yfun(Yn_1, b2, betan, delSn, A2, delalphan) #x
    return Yn_1*exp(-b2*(betan^2)*delSn) + A2*delalphan*exp(-b2*(betan^2)*delSn/2)
end

function Xfun(Xn_1, b1, betan, delSn, A1, delalphan) #x
    return Xn_1*exp(-b1*(betan^2)*delSn) + A1*delalphan*exp(-b1*(betan^2)*delSn/2)
end

function alphaEfun(alphan, Xn, Yn) #x
    return alphan-Xn-Yn
end

function CcNfun(CNalpha, alphaEn) #x
    return CNalpha*alphaEn
end

function kalphafun(Mn, betan, A1, b1, A2, b2) #x #Associated with the noncirculatory time constant
    return 0.75/((1-Mn) + (pi*(betan^2)*(Mn^2)*(A1*b1 + A2*b2)))
end

function Dfun(Dn_1, delt, kalphan, Tl, delalphan, delalphan_1) #x
    t1 = Dn_1*exp(-delt/(kalphan*Tl))
    t2 = ((delalphan-delalphan_1)/delt)*exp(-delt/(2*kalphan*Tl))
    return t1+t2
end

function CINfun(kalphan, Tl, Mn, delalphan, delt, Dn) #x
    return 4*kalphan*Tl*((delalphan/delt)-Dn)/Mn
end

function CPNfun(CINn, CCNn) #x
    return CINn + CCNn
end

function DPfun(DPn_1, delS, Tp, CPNn, CPNn_1) #x
    t1 = DPn_1*exp(delS/Tp)
    t2 = (CPNn-CPNn_1)*exp(delS/(2*Tp))
    return t1 + t2
end

#Apparently Tp is a function of M

function CNfun(alphan, alpha0, Cnalpha, f) #x
    return Cnalpha*(((1 + sqrt(f))/2)^2)*(alphan-alpha0) #Note: This is different then the original formulation. In the original they multiply by alphan. 
end

function CprimeNfun(CNn, DPn) #x
    return CNn-DPn
end

function alphaffun(CprimeN, CNalpha) #x
    return CprimeN/CNalpha
end

function fprimefun(alphafn, alpha1, S1, S2) #x
    return ffun(alphafn, alpha1, S1, S2)
end

function fppfun(fprimen, Dfn) #x
    return fprimen-Dfn
end

function CfNfun(CNalpha, fppn, alphaEn, CINn) #x
    # println("fpp factor: ", ((((1+sqrt(fppn))/2))^2))
    return CNalpha*((((1+sqrt(fppn))/2))^2)*alphaEn + CINn
end


function BeddoesLeishman(tspan, alpha, U, c, A, b, S, alpha0, alpha1, CNalpha;Tp=1.7, Tf=3.0, Tv=6.0, nt=100, f0=NaN, a=343.0)
    # tspan - range of time to calculate across
    # alpha - angle of attack function across time
    # U - X direction windspeed function across time
    # c - chord value
    # A - 
    # b - 
    # S - Static stall characteristics - "Easily determined from the static lift data" - pg 7
    # alpha0 - Zero lift angle of attack
    # alpha1 - Static Stall angle of attack
    # CNalpha - attached (linear) region slope - Can put as 2*pi
    # nt - number of time steps
    # f0 - Initial separation value
    # a - speed of sound
    


    ### Find the time spacing
    tvec = collect(range(tspan[1], tspan[2]; length=nt))

    ### Environmental Constants
    deltat = tvec[2]-tvec[1] #Assuming that the second time step is the same length as the first. 
    M = zeros(nt)
    beta = zeros(nt)

    ### Angle of attack values
    alphaf = zeros(nt)
    alphaE = zeros(nt)
    delalpha = zeros(nt)
    #Need alpha1

    ### Separation Factor Values
    f = zeros(nt)
    fp = zeros(nt)
    fpp = zeros(nt)

    ### Deficiency Function values
    D = zeros(nt)
    DP = zeros(nt)
    Df = zeros(nt)
    X = zeros(nt)
    Y = zeros(nt)

    ### Intermediate Normal force values (Normal Force is Lift)
    CfN = zeros(nt)
    CvN = zeros(nt)
    CN = zeros(nt)
    CPN = zeros(nt)
    CprimeN = zeros(nt)
    CIN = zeros(nt)
    CcN = zeros(nt)
    Cv = zeros(nt)
    CvN = zeros(nt)

    ### Outputs
    Clvec = zeros(nt)
    Cdvec = zeros(nt)
    CNdyn = zeros(nt)

    ### Calculated Constants
    kalpha = zeros(nt)
    delS = zeros(nt)
    KN = zeros(nt)

    ### Airfoil Constants
    #Need CNalpha
    Tl = c/a
    #Tf
    #Tp
    #Tv
    #Tvl? - I don't use this, but it is in the paper.... meep. 
    ## Name year (Tp, Tf, Tv, Tvl)
    # Leishman 86 (1.7, 3.0, 6.0, 7.0)
    # Leishman 89 (2.0, 2.5, 6.0, 11.0)
    # Pierce 96 (1.7, 3.0, 6.0, 11.0)
    # Minnema 98 (1.7, 3.0, 6.0, 11.0)
    ## Coefficients of separation point curve fit -> I think they are just coefficients to make the inviscid lift match the static lift (probably in the attached region). 
    S1 = S[1] #Wind Energy Handbook gives a value of 0.05 for S1 and S2
    S2 = S[2]
    A1 = A[1]
    A2 = A[2]
    b1 = b[1]
    b2 = b[2]

    ### Initial Values
    M[1] = U(tvec[1])/a
    beta[1] = betafun(M[1])

    if isnan(f0)
        f[1] = ffun(alpha(tvec[1]), alpha1, S1, S2)
    else
        f[1] = f0
    end

    CN[1] = CNfun(alpha(tvec[1]), alpha0, CNalpha, f[1]) #x
    CprimeN[1] = CprimeNfun(CN[1], DP[1]) #x
    alphaf[1] = alphaffun(CprimeN[1], CNalpha) #x
    fp[1] = fprimefun(alphaf[1], alpha1, S1, S2) #x
    fpp[1] = fppfun(fp[1], Df[1]) #x

    alphaE[1] = alphaEfun(alpha(tvec[1]), X[1], Y[1]) #x
    kalpha[1] = kalphafun(M[1], beta[1], A1, b1, A2, b2) #x
    CIN[1] = CINfun(kalpha[1], Tl, M[1], delalpha[1], deltat, D[1]) #x #Implusive lift
    CcN[1] = CcNfun(CNalpha, alphaE[1]) #x
    CPN[1] = CIN[1] + CcN[1] #x

    CcN[1] = CcNfun(CNalpha, alphaE[1]) #x
    KN[1] = KNfun(fpp[1]) #x
    Cv[1] = Cvfun(KN[1], CcN[1]) #x #vortex lift 

    CfN[1] = CfNfun(CNalpha, fpp[1], alphaE[1], CIN[1])
    CvN[1] = CvNfun(CvN[1], Cv[1], Cv[1], delS[1], Tv)


    CNdyn[1] = CfN[1] + CvN[1]


    #Todo: Make sure all the constants are calculated/retrieved
    for i = 2:nt
        deltat = tvec[i]-tvec[i-1]
        alphan = alpha(tvec[i]) #x
        delalpha[i] = alphan-alpha(tvec[i-1]) #x - Letting it remain zero
        delS[i] = delSfun(c, U(tvec[i]), U(tvec[i-1]), tvec[i], tvec[i-1]) #x - letting it remain zero # S - non-dimensional distance traveled
        M[i] = U(tvec[i])/a #x
        beta[i] = betafun(M[i]) #x

        ### Find CfN
        X[i] = Xfun(X[i-1], b1, beta[i], delS[i], A1, delalpha[i]) #x - i.v.=0
        Y[i] = Yfun(Y[i-1], b2, beta[i], delS[i], A2, delalpha[i]) #x - i.v.=0
        kalpha[i] = kalphafun(M[i], beta[i], A1, b1, A2, b2) #x - i.v. unused
        D[i] = Dfun(D[i-1], deltat, kalpha[i], Tl, delalpha[i], delalpha[i-1]) #x - i.v.=0
        f[i] = ffun(alphan, alpha1, S1, S2) #x
        CN[i] = CNfun(alphan, alpha0, CNalpha, f[i]) #x
        CIN[i] = CINfun(kalpha[i], Tl, M[i], delalpha[i], deltat, D[i]) #x #Implusive lift
        alphaE[i] = alphaEfun(alphan, X[i], Y[i]) #x #initial values for X, Y are 0.
        CcN[i] = CcNfun(CNalpha, alphaE[i]) #x
        CPN[i] = CIN[i] + CcN[i] #x
        DP[i] = DPfun(DP[i-1], delS[i], Tp, CPN[i], CPN[i-1]) # Todo: Not sure what to do about Tp. One source gives a value for it. The original paper says that it is a function of Mach number, but I can't find the equation or an approximation. 
        CprimeN[i] = CprimeNfun(CN[i], DP[i]) #x
        alphaf[i] = alphaffun(CprimeN[i], CNalpha) #x
        fp[i] = fprimefun(alphaf[i], alpha1, S1, S2) #x
        Df[i] = Dffun(Df[i-1], delS[i], Tf, fp[i], fp[i-1]) #Todo: similar problem with Tf. However there is a source that claims that you can caluclate it. 
        fpp[i] = fppfun(fp[i], Df[i]) #x
        CfN[i] = CfNfun(CNalpha, fpp[i], alphaE[i], CIN[i]) #x

        ### Find CvN
        KN[i] = KNfun(fpp[i]) #x
        Cv[i] = Cvfun(KN[i], CcN[i]) #x #vortex lift 
        CvN[i] = CvNfun(CvN[i-1], Cv[i], Cv[i-1], delS[i], Tv) #Todo: Tv? others done.  #Total Vortex lift - Assuming the initial condition is attached, thus vortex lift should be zero. 
        
        #Todo: This is for dynamic stall conditions, what about in steady flow? 
        println("$i th iteration")
        println("beta: ", beta[i])
        println("alpha: ", alphan)
        println("delS: ", delS[i])
        println("X[i]: ", X[i])
        println("Y[i]: ", Y[i])
        println("D: ", D[i])
        println("DP: ", DP[i])
        println("Df: ", Df[i])
        println("f: ", f[i])
        println("fp: ", fp[i])
        println("fpp: ", fpp[i])
        println("alphaE: ", alphaE[i])
        println("CN: ", CN[i])
        println("CIN: ", CIN[i])
        println("CcN: ", CcN[i])
        println("Cpn: ", CPN[i])
        println("Cfn : ", CfN[i])
        println("CvN : ", CvN[i])
        println("")

        # if X[i]>0
        #     println("X[$i]: ", X[i])
        # end
        CNdyn[i] = CfN[i] + CvN[i]
        
    end
    return CNdyn, tvec
end