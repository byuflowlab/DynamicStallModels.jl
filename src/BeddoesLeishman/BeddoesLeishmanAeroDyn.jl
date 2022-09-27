#=
AeroDyn's implementation of the Beddoes-Leishman model. 

=#


function getloads_BLA(dsmodel::BeddoesLeishman, states, p, airfoil)
    c, a, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, alpha1, alpha2, S1, S2, S3, S4, xcp, U, _ = p
    Cnfit = airfoil.cl

    Cn, Cc, Cl, Cd, Cm = BLAD_coefficients(dsmodel::BeddoesLeishman, states, U, c, Cnfit, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tvl, xcp, a)
    return [Cn, Cc, Cl, Cd, Cm]
end



function seperationpoint(alpha, alpha0, alpha1, alpha2, S1, S2, S3, S4)
    #Question: I think alpha1==afp and alpha2==afm
    if (alpha0 <= alpha < alpha1)
        return 1 - 0.3*exp((alpha-alpha1)/S1)
    elseif (alpha2 <= alpha < alpha0)
        return 1 - 0.3*exp((alpha2-alpha)/S3)
    elseif alpha>alpha1
        return 0.04 + 0.66*exp((alpha1-alpha)/S2)
    elseif alpha<alpha2
        return 0.04 + 0.66*exp((alpha-alpha2)/S4)
    end

    @warn("No seperation point found. Returning 1.0.")
    return 1.0
end


#AeroDyn original implementation. 
function update_states_ADO(dsmodel::BeddoesLeishman, oldstates, c, a, U, deltat, aoa, dcndalpha, alpha0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, alpha1, alpha2, S1, S2, S3, S4) 

    ### Unpack
    alpha_m, alphaf_m, q_m, Ka_m, Kq_m, X1_m, X2_m, Kpa_m, Kpq_m, Kppq_m, Kpppq_m, Dp_m, Df_m, Cpotn_m, fp_m, fpp_m, tauv, Cvn_m, Cv_m, LESF_m, TESF_m, VRTX_m = oldstates #The underscore m means that it is the previous time step (m comes before n).

    # if c==4.557
    #     @show aoa
    # end

    #= States
    1 - alpha
    2 - alphaf #TODO: I'm not sure that this state is needed. 
    3 - q
    4 - Ka
    5 - Kq
    6 - X1
    7 - X2
    8 - Kpa
    9 - Kpq
    10 - Kppq
    11 - Kpppq
    12 - Dp
    13 - Df
    14 - Cpotn
    15 - fp
    16 - fpp
    17 - tauv
    18 - Cvn
    19 - Cv
    20 - LESF::Bool
    21 - TESF::Bool
    22 - VRTX::Bool
    =#


    states = zeros(22)




    zeta, A5, b5, Tsh, _ = dsmodel.constants 
    #=
    zeta - Low-pass-fileter frequency cutoff. #TODO: Should this have the negative or should the equation have the negative? 
    Tsh - Strouhal's frequency 0.19
    =#

    ########### Algorithm ############### (Converted from UA documentation)
    ### Initial constants
    TI = c/a # Equation 1.11
    deltas = 2*U*deltat/c# Equation 1.5b



    ###Low-pass-filtering #TODO: Might put this in a function. lowpass()
    Clp = exp(2*pi*deltat*zeta) #low-pass-filter constant. #Equation 1.8g -> removed the negative because zeta has a negative. 
    plC = 1 - Clp

    states[1] = alpha = Clp*alpha_m + plC*aoa #low-pass-filtered angle of attack. #Equation 1.8a
    # @show alpha

    q = (alpha - alpha_m)*c/(U*deltat) #Pitch rate. #Equation 1.8b #This appears to be dimensionless. 
    states[3] = q = Clp*q_m + plC*q #low-pass-filter pitch rate. #Equation 1.8c

    states[4] = Ka = q*U/deltat #Equation 1.8d #Todo: this combination of q and delta t don't work for small delta t. It makes Ka super larger. So either q needs to be smaller, or delta t larger. If delta t is larger, then q should be smaller as well... which is weird, because that puts a bottom limit on delta t... which breaks CFD rules. What is it Patankar's rules? I think this breaks Patankar's rules. 
    # @show q, deltat #because deltat is small, Ka is large.... Maybe q should be smaller? 

    # if c== 1.419
    #     @show q, U, deltat
    # end

    Kq = (q - q_m)/deltat #Equation 1.8e
    states[5] = Kq = Clp*Kq_m + plC*Kq #Equation 1.8f

    ### Update sigma 1 - Tf modifications
    Delta_alpha0 = alpha - alpha0
    sigma1 = 1
    sigma3 = 1
    if TESF_m == 1 #(Separation)
        if Ka_m*Delta_alpha0 < 0
            sigma1 = 2 #(Accelerate separation point movement)
        else
            if LESF_m == 0
                sigma1 = 1 #(LE separation can occur)
            else
                if fpp_m <= 0.7 
                    sigma1 = 2 #(accelerate separation point movement if separation is occuring )
                else
                    sigma1 = 1.75
                end
            end
        end
    else #(reattachment (`TESF==false`))
        if LESF_m == 0
            sigma1 = 0.5 #(Slow down reattachment)
        elseif VRTX_m == 1  && 0 <= tauv <= Tvl
            sigma1 = 0.25 #(No flow reattachment if vortex shedding is in progress)
        elseif Ka_m*Delta_alpha0 > 0
            sigma1 = 0.75
        end 
    end

    ### Update sigma3 - Tv modifications
    if Tvl <= tauv <= 2*Tvl
        sigma3 = 3 #Postshedding
        if TESF_m == 0
            sigma3 = 4 #Accelerate vortex lift decay
            if VRTX_m == 1 && 0 <= tauv <= Tvl
                if Ka_m*Delta_alpha0 < 0
                    sigma3 = 2 #Accelerate vortex lift decay
                else
                    sigma3 = 1 #default
                end
            end
        end
    else
        if Ka_m*Delta_alpha0 < 0
            sigma3 = 4 #vortex lift must decay fast
        end
    end
    
    if TESF_m == 0 && Kq_m*Delta_alpha0 < 0 #Note that it's Kq and not K_alpha
        sigma3 = 1 #Default
    end



    ### More constants
    M = U/a # Mach number
    beta = sqrt(1 - M^2) #Prandtl-Glauert compressibility correction factor #TODO: Need something to cap if M>1. Or at least some sort of diversion behavior. 

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) #TODO: Calculated solely at first time step.... Maybe we pull this out and pass it in as a function argument. ... Or just calculate it every time as it isn't super costly. 
    k_alpha = 1/bot #Equation 1.11a

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) #TODO. The documentation has the same equation for 1.11b as 1.11a. -> In Leishman's 1900 state space paper, he has the second term in the denominator of k_alpha to be half of what is shown here. I'll try what I have, and if it is off, then I'll try changing it. -> It appears to be correct. 
    k_q = 1/bot #Equation 1.11b

    Talpha = 3*k_alpha*TI/4 #Equation 1.10a
    Tq = 3*k_q*TI/4 #Equation 1.10b

    Tf = Tf0/sigma1 #Equation 1.37
    Tv = Tv0/sigma3 #Equation 1.48



    ### Non-circulatory components of normal force
    states[8] = Kpa = Kpa_m*exp(-deltat/Talpha) + (Ka - Ka_m)*exp(-deltat/(2*Talpha)) #Equation 1.18b, Deficiency function for eq. 1.18 
    Cnc_nalpha = 4*Talpha*(Ka-Kpa)/M #Equation 1.18, Noncirculatory component of normal force due to changes in alpha

    states[9] = Kpq = Kpq_m*exp(-deltat/Tq) + (Kq - Kq_m)*exp(-deltat/(2*Tq)) #Equation 1.19b, deficiency function for eq. 1.19
    Cnc_nq = Tq*(Kq - Kpq)/M #Equation 1.19, Noncirculatory component of normal force due to changes in pitching rate. #TODO. The documentation conflicts on whether or not a minus should be included here. -> I'm going to assume that it is positive. And if I'm wrong... I'll change it. It seems to be correct with positive. 

    Cnc_naq = Cnc_nalpha + Cnc_nq #Equation 1.17, Noncirculatory component of normal force via superposition. 


    




    ### Update states 1 and 2
    beta2 = beta^2 #beta squared
    delta_alpha = alpha - alpha_m #Note: Not explicitly defined in the documentation. Assumed. 

    states[6] = X1 = X1_m*exp(-b1*beta2*deltas) + A1*exp(-b1*beta2*deltas/2)*delta_alpha #EQ 1.15a
    states[7] = X2 = X2_m*exp(-b2*beta2*deltas) + A2*exp(-b2*beta2*deltas/2)*delta_alpha #EQ 1.15b



    ### Circulatory component of normal force
    alphae = (alpha - alpha0) - X1 - X2 #EQ 1.14, Effective angle of attack

    Cc_na = dcndalpha/beta #EQ 1.12, Circulatory component of the normal force coefficient response to step change in alpha. 
    Cc_naq = Cc_na*alphae #EQ 1.13, Circulatory component of normal force via lumped approach. #TODO: This appears to be equal to Cpotcn



    ### Circulatory component of moment. #Question: Why did they calculate the moment here? Do I need some of these things here? Is there a better spot to put this? 
    deltaq = q - q_m #Not explicityly stated in the docs. Assumed
    states[11] = Kpppq = Kpppq_m*exp(-b5*beta2*deltas) + A5*deltaq*exp(-b5*beta2*deltas/2) #EQ 1.26



    ### Total normal force under attached conditions
    states[14] = Cpotn = Cc_naq + Cnc_naq #Equation 1.20
    # if c== 1.419
    #     @show Cc_naq, Cnc_naq
    # end


    ### Noncirculatory component of moment due to change in pitch #Note: If this can, this should probably go with the other moment calculation. 
    bot = 15*(1-M) + 3*dcndalpha*A5*b5*beta*M*M/2  
    kmq = 7/bot #EQ 1.29b
    states[10] = Kppq = Kppq_m*exp(-deltat/(kmq*kmq*TI)) + (Kq - Kq_m)*exp(-deltat/(2*kmq*kmq*TI)) #EQ 1.29c



    ####### boundary layer response 
    states[12] = Dp = Dp_m*exp(-deltas/Tp) + (Cpotn - Cpotn_m)*exp(-deltas/(Tp*2)) #EQ1.35b, deficiency function. 
    Cpn = Cpotn - Dp #EQ 1.35, lagged circulatory normal force 


    states[2] = alphaf = Cpn/Cc_na + alpha0 #EQ 1.34, delayed effective angle of incidence  

    states[15] = fp = seperationpoint(alphaf, alpha0, alpha1, alpha2, S1, S2, S3, S4) #EQ 1.33 
    states[13] = Df = Df_m*exp(-deltas/Tf) + (fp - fp_m)*exp(-deltas/(2*Tf)) #EQ 1.36b
    states[16] = fpp = fp - Df #EQ 1.36, delayed effective seperation point. 
    


    ### 
    fterm = (1 + sqrt(fpp))/2



    states[19] = Cv = Cc_na*alphae*(1 - fterm)^2 #EQ 1.49, Normal force coefficient due to accumulated vorticity

    states[18] = Cvn = Cvn_m*exp(-deltas/Tv) + (Cv - Cv_m)*exp(-deltas/(2*Tv)) #EQ 1.47
    #Note: Cv has to be the same sign as Cfs_n
    # if c==1.419
    #     @show Cvn, Cc_na, alphae, fterm
    # end
    


    ######## Update "other states"
    ### Test for Leading edge separation
    if Cpn > Cn1 
        LESF = 1.0 #LE separation can occur
    else
        LESF = 0.0 #Reattachment can occur
    end

    ### Test for Trailing edge separation
    if fpp < fpp_m
        TESF = 1.0 #TE separation in progress
    else
        TESF = 0.0
    end

    ### Test for vortex advection
    if 0< tauv <= 2*Tvl
        VRTX = 1.0 #Vortex advection in progress
    else
        VRTX = 0.0 #Vortex is in wake. 
    end

    ### Vortex position reset
    if (tauv >= 1 + Tsh/Tvl) & (LESF == 1.0)
        tauv = 0
    end

    if LESF==1
        tauv += deltat*2*U/c #No way given to update tauv. Doing a simple Euler step. 
    end

    

    states[17] = tauv
    states[20] = LESF
    states[21] = TESF
    states[22] = VRTX



    return states
end

function BLAD_coefficients(dsmodel::BeddoesLeishman, states, U, c, af::Airfoil, a)
    cnfit = af.cl
    dcndalpha = af.dcldalpha
    alpha0 = af.alpha0
    Cd0 = af.cd(alpha0)
    Cm0 = af.cm(alpha0)
    A1 = af.A[1]
    A2 = af.A[2]
    b1 = af.b[1]
    b2 = af.b[2]
    Tvl = af.T[4]
    xcp = af.xcp

    return BLAD_coefficients(dsmodel, states, U, c, cnfit, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tvl, xcp, a)
end

function BLAD_coefficients(dsmodel::BeddoesLeishman, states, U, c, Cnfit, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tvl, xcp, a)

    _, A5, b5, _, eta = dsmodel.constants 

    Ka = states[4]
    Kpa = states[8]
    Kq = states[5]
    Kpq = states[9]

    M = U/a # Mach number
    beta = sqrt(1 - M^2) #Prandtl-Glauert compressibility correction factor
    TI = c/a




    ##### Prepare inputs for normal force coefficient. 
    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2)  
    k_alpha = 1/bot #Equation 1.11a
    Talpha = 3*k_alpha*TI/4 #Equation 1.10a
    Cnc_nalpha = 4*Talpha*(Ka-Kpa)/M #Equation 1.18, Noncirculatory component of normal force due to changes in alpha

    # if c== 1.419
    #     @show Talpha, Ka, Kpa, M #Ka and Kpa are quite large. 
    # end

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) 
    k_q = 1/bot #Equation 1.11b
    Tq = 3*k_q*TI/4 #Equation 1.10b
    Cnc_nq = Tq*(Kq - Kpq)/M #Equation 1.19

    Cnc_naq = Cnc_nalpha + Cnc_nq #Equation 1.17, Noncirculatory component of normal force via superposition.

    # if c== 1.419
    #     @show Cnc_nalpha, Cnc_nq #Cnc_nalpha is quite large. 
    # end

    alpha = states[1]
    X1 = states[6]
    X2 = states[7]
    alphae = (alpha - alpha0) - X1 - X2 #EQ 1.14, Effective angle of attack
    Cc_na = dcndalpha/beta #EQ 1.12, Circulatory component of the normal force coefficient response to step change in alpha. 
    Cc_naq = Cc_na*alphae #EQ 1.13, Circulatory component of normal force via lumped approach.

    fpp = states[16]
    fterm = (1 + sqrt(fpp))/2

    Cfsn = Cnc_naq + Cc_naq*(fterm)^2 #EQ 1.38, Normal force coefficient after accounting for separated flow from TE
    # if c== 1.419
    #     @show Cnc_naq, Cc_naq, fterm #Both the Cnc and the Cc terms are large, the Cnc term is especially large. 
    # end

    Cvn = states[18]

    ### Total normal force 
    Cn = Cfsn + Cvn #EQ 1.53

    # if c== 1.419
    #     @show Cfsn, Cvn #Problem is Cfsn
    # end



    ######### Prepare inputs for chordwise force coefficient. 
    Cpotc = Cc_naq*tan(alphae + alpha0) #Equation 1.21  
    Cfsc = Cpotc*eta*sqrt(fpp) #EQ 1.40 #(sqrt(fpp)-0.2) #Gonzalez modifications

    ### Chordwise force 
    Cc = Cfsc + Cvn*tan(alphae) #EQ 1.55




    ### Lift and Drag
    Cl = Cn*cos(alpha) + Cc*sin(alpha)
    Cd = Cn*sin(alpha) - Cc*cos(alpha) + Cd0





    q = states[3]
    Kpppq = states[11]
    Cc_mq = -dcndalpha*(q-Kpppq)*c/(16*beta*U) #Equation 1.22c

    ### Noncirculatory component of moment due to change in alpha
    Cnc_ma = -Cnfit(alpha)/4 #EQ 1.27 Note: This equation was missed in the documented algorithm.

    bot = 15*(1-M) + 3*dcndalpha*A5*b5*beta*M*M/2  
    kmq = 7/bot #EQ 1.29b
    Kppq = states[10]
    Cnc_mq = -7*TI*kmq*kmq*(Kq-Kppq)/(12*M) #EQ 1.29, Circulatory component of moment.

    tauv = states[17]
    xbarcp = 0.2 #Todo: Is this different from xcp? They call it x bar bar cp
    xvcp = xbarcp*(1-cos(pi*tauv/(Tvl))) #1.57b
    Cvm = -xvcp*Cvn #1.57a

    ### Moment
    Cm = Cm0 - Cc_naq*(xcp - 0.25) + Cc_mq + Cnc_ma + Cnc_mq + Cvm #Equation 1.58


    return Cn, Cc, Cl, Cd, Cm
end

function initialize_ADO(Uvec, aoavec, tvec, airfoil::Airfoil, c, a) 
    # Cnfit = airfoil.cl
    dcndalpha = airfoil.dcldalpha
    alpha0 = airfoil.alpha0
    A1, A2 = airfoil.A
    b1, b2 = airfoil.b
    Tp, Tf0, Tv0, Tvl = airfoil.T
    alpha2, alpha1 = airfoil.alphasep
    S1, S2, S3, S4 = airfoil.S
    xcp = airfoil.xcp

    Cn1 = airfoil.cl(alpha1)
    Cd0 = airfoil.cd(alpha0)
    Cm0 = airfoil.cm(alpha0)



    dt = tvec[2] - tvec[1]

    aoa = aoavec[1]
    alpha = alphaf = aoavec[1]
    aoadot = q = Ka = 0 #(aoavec[2] - aoavec[1])/dt #Ka was initialized too high. aoadot isn't used. I might be able to use q = (aoavec[2] - aoavec[1])*c/(U*deltat) and ka = q*U/deltat... Why are q and Ka both states?... if one is just a multiple of the other. I guess U and delta t change, so it isn't the same multiple across time. 
    Kq = 0.0
    X1 = X2 = 0.0
    Kpa = Kpq = Kppq = Kpppq = 0.0
    Dp = Df = 0.0
    Cpotn = airfoil.cl(aoavec[1])
    fp = fpp = 1.0
    tauv = 0.0
    Cvn = 0.0
    Cv = 0.0

    LESF = TESF = VRTX = 0

    states = [alpha, alphaf, q, Ka, Kq, X1, X2, Kpa, Kpq, Kppq, Kpppq, Dp, Df, Cpotn, fp, fpp, tauv, Cvn, Cv, LESF, TESF, VRTX]

    loads = [airfoil.cl(aoa), airfoil.cd(aoa), airfoil.cl(aoa), airfoil.cd(aoa), airfoil.cm(aoa)]


    p = [c, a, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, alpha1, alpha2, S1, S2, S3, S4, xcp]

    envvars = [Uvec[1], aoavec[1]]

    return states, loads, vcat(p, envvars)
end

function updateenvironment_ADO(p, U, aoa)
    p[23] = U
    p[24] = aoa
end