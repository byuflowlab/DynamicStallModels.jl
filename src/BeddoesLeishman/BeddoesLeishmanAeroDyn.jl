#=
AeroDyn's implementation of the Beddoes-Leishman model. 

=#

function (model::BeddoesLeishman)(x, p, y)
    if isa(model.detype, Functional)
        @warn("Functional implementation not yet prepared.")
    elseif isa(model.detype, Indicial)
        if model.version==1
            @warn("Original indicial Beddoe-Leishman not prepared for use yet.")
        elseif model.version==2
    
            c, a, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, alpha1, alpha2, S1, S2, S3, S4, xcp, _, _, _ = p #Inputs
            U, aoa, dt = y #Environmental inputs. 

            flags = view(p, 23:25)

            return update_states_ADO(model, x, flags, c, a, U, dt, aoa, dcndalpha, alpha0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, alpha1, alpha2, S1, S2, S3, S4)

        elseif model.version==3
            @warn("AeroDyn Beddoe-Leishman with Gonzalez's modifications not prepared for use yet.")
        elseif model.version==4
            @warn("AeroDyn Beddoe-Leishman with Minema's modifications not prepared for use yet.")
        end
    end
    
end

function getloads(dsmodel::BeddoesLeishman, states, p, y, airfoil)
    c, a, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, alpha1, alpha2, S1, S2, S3, S4, xcp, _, _, _ = p
    Cnfit = airfoil.cl

    U, _, _ = y

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

    @warn("No seperation point found. Return 1.0.")
    return 1.0
end


#AeroDyn original implementation. 
function update_states_ADO(dsmodel::BeddoesLeishman, oldstates, flags, c, a, U, deltat, aoa, dcndalpha, alpha0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, alpha1, alpha2, S1, S2, S3, S4)

    ### Unpack
    aoa_m, alpha_m, alphaf_m, aoadot_m, q_m, Ka_m, Kq_m, X1_m, X2_m, X3_m, X4_m, Kpa_m, Kpq_m, Kppq_m, Kpppq_m, Dp_m, Df_m, Dfc_m, Cpotn_m, fp_m, fpc_m, fpp_m, fppc_m, tauv, Cvn_m, Cv_m, Daf_m, sigma1, sigma3 = oldstates #The underscore m means that it is the previous time step (m comes before n).

    LESF, TESF, VRTX = flags

    states = zeros(29)

    
    states[1] = aoa
    states[4] = aoadot = 0 #Appears to be unused. 
    states[10] = X3 = 0 #Appears to be unused.
    states[11] = X4 = 0 #Appears to be unused.
    states[18] = Dfc = 0 #Appears to be unused.
    states[21] = fpc = 1 #Appears to be unused.
    states[23] = fppc = 1 #Appears to be unused.
    # states[24] = tauv = 0 #Appears to be unused.
    states[27] = Daf = 0 #Appears to be unused.

    #Removal of states? 
    #=
    - aoa_m - The angle of attack at the previous time step (before low-pass filtering). )
    - alphaf_m - delayed effective angle of incidence. 
    - aoadot_m
    - X3
    - X4
    - Dfc
    - fpc 

    =#


    zeta, A5, b5, Tsh, _ = dsmodel.constants 
    #=
    zeta - Low-pass-fileter frequency cutoff. #TODO: Should this have the negative or should the equation have the negative? 
    Tsh - Strouhal's frequency 0.19

    
    =#

    ########### Algorithm ############### (Converted from UA documentation)
    ### Initial constants
    TI = c/a # Equation 1.11
    deltas = 2*U*deltat/c# Equation 1.5b



    # Kan = (aoa - aoa_m)/deltat
    # q = Kan*c/U # Equation 1.7 #Question: I don't think that equation 1.7 ever actually gets used. -> Maybe it's supposed to be an if statement? I dunno. 



    ###Low-pass-filtering #TODO: Might put this in a function. lowpass()
    Clp = exp(2*pi*deltat*zeta) #low-pass-filter constant. #Equation 1.8g -> removed the negative because zeta has a negative. 
    plC = 1 - Clp

    states[2] = alpha = Clp*alpha_m + plC*aoa #low-pass-filtered angle of attack. #Equation 1.8a
    # @show alpha

    q = (alpha - alpha_m)*c/(U*deltat) #Pitch rate. #Equation 1.8b #This appears to be dimensionless. 
    states[5] = q = Clp*q_m + plC*q #low-pass-filter pitch rate. #Equation 1.8c

    states[6] = Ka = q*U/deltat #Equation 1.8d #Todo: this combination of q and delta t don't work for small delta t. It makes Ka super larger. So either q needs to be smaller, or delta t larger. If delta t is larger, then q should be smaller as well... which is weird, because that puts a bottom limit on delta t... which breaks CFD rules. What is it Patankar's rules? I think this breaks Patankar's rules. 
    # @show q, deltat #because deltat is small, Ka is large.... Maybe q should be smaller? 

    Kq = (q - q_m)/deltat #Equation 1.8e
    states[7] = Kq = Clp*Kq_m + plC*Kq #Equation 1.8f



    ### More constants
    M = U/a # Mach number
    beta = sqrt(1 - M^2) #Prandtl-Glauert compressibility correction factor #TODO: Need something to cap if M>1. Or at least some sort of diversion behavior. 

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) #TODO: Calculated solely at first time step.... Maybe we pull this out and pass it in as a function argument. ... Or just calculate it every time as it isn't super costly. 
    k_alpha = 1/bot #Equation 1.11a

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) #TODO: The documentation has the same equation for 1.11b as 1.11a. -> In Leishman's 1900 state space paper, he has the second term in the denominator of k_alpha to be half of what is shown here. I'll try what I have, and if it is off, then I'll try changing it. 
    k_q = 1/bot #Equation 1.11b

    Talpha = 3*k_alpha*TI/4 #Equation 1.10a
    Tq = 3*k_q*TI/4 #Equation 1.10b

    Tf = Tf0/sigma1 #Equation 1.37
    Tv = Tv0/sigma3 #Equation 1.48



    ### Non-circulatory components of normal force
    states[12] = Kpa = Kpa_m*exp(-deltat/Talpha) + (Ka - Ka_m)*exp(-deltat/(2*Talpha)) #Equation 1.18b, Deficiency function for eq. 1.18 
    Cnc_nalpha = 4*Talpha*(Ka-Kpa)/M #Equation 1.18, Noncirculatory component of normal force due to changes in alpha
    # @show Talpha, Ka, Kpa, M #Both Ka and Kpa are large. I think Kpa is large because Ka is large. 

    states[13] = Kpq = Kpq_m*exp(-deltat/Tq) + (Kq - Kq_m)*exp(-deltat/(2*Tq)) #Equation 1.19b, deficiency function for eq. 1.19
    Cnc_nq = Tq*(Kq - Kpq)/M #Equation 1.19, Noncirculatory component of normal force due to changes in pitching rate. #TODO: The documentation conflicts on whether or not a minus should be included here. -> I'm going to assume that it is positive. And if I'm wrong... I'll change it. 

    Cnc_naq = Cnc_nalpha + Cnc_nq #Equation 1.17, Noncirculatory component of normal force via superposition. 

    # @show Cnc_nalpha, Cnc_nq #Cnc_nalpha is super off. Again. 

    




    ### Update states 1 and 2
    beta2 = beta^2 #beta squared
    delta_alpha = alpha - alpha_m #Note: Not explicitly defined in the documentation. Assumed. 

    states[8] = X1 = X1_m*exp(-b1*beta2*deltas) + A1*exp(-b1*beta2*deltas/2)*delta_alpha #EQ 1.15a
    states[9] = X2 = X2_m*exp(-b2*beta2*deltas) + A2*exp(-b2*beta2*deltas/2)*delta_alpha #EQ 1.15b



    ### Circulatory component of normal force
    alphae = (alpha - alpha0) - X1 - X2 #EQ 1.14, Effective angle of attack

    Cc_na = dcndalpha/beta #EQ 1.12, Circulatory component of the normal force coefficient response to step change in alpha. 
    Cc_naq = Cc_na*alphae #EQ 1.13, Circulatory component of normal force via lumped approach. #TODO: This appears to be equal to Cpotcn



    ### Circulatory component of moment. #Question: Why did they calculate the moment here? Do I need some of these things here? Is there a better spot to put this? 
    deltaq = q - q_m #Not explicityly stated in the docs. Assumed
    states[15] = Kpppq = Kpppq_m*exp(-b5*beta2*deltas) + A5*deltaq*exp(-b5*beta2*deltas/2) #EQ 1.26
    # Cc_mq = -dcndalpha*(q-Kpppq)*c/(16*beta*U) #EQ 1.25, Circulatory component of moment


    ### Total normal force under attached conditions
    states[19] = Cpotn = Cc_naq + Cnc_naq #Equation 1.20
    # @show Cc_naq, Cnc_naq #Cnc_naq is gigantic. 


    ### Noncirculatory component of moment due to change in alpha
    # Cnc_ma = -Cnfit(alpha)/4 #EQ 1.27 Note: This equation was missed in the documented algorithm. 



    ### Noncirculatory component of moment due to change in pitch #Note: If this can, this should probably go with the other moment calculation. 
    bot = 15*(1-M) + 3*dcndalpha*A5*b5*beta*M*M/2  
    kmq = 7/bot #EQ 1.29b
    states[14] = Kppq = Kppq_m*exp(-deltat/(kmq*kmq*TI)) + (Kq - Kq_m)*exp(-deltat/(2*kmq*kmq*TI)) #EQ 1.29c
    # Cnc_mq = -7*TI*kmq*kmq*(Kq-Kppq)/(12*M) #EQ 1.29




    ### Chordwise force
    # Cpotcn = Cc_na*alphae #Todo: I think that this C^c_{n\alpha}(s,M) should be a circulatory normal coefficient. ... Well.. Now I don't know. It could be the compressibility corrected slope. 
    # Cpotc = Cpotcn*tan(alphae + alpha0) #Equation 1.21  



    ####### boundary layer response 
    states[16] = Dp = Dp_m*exp(-deltas/Tp) + (Cpotn - Cpotn_m)*exp(-deltas/(Tp*2)) #EQ1.35b, deficiency function. 
    Cpn = Cpotn - Dp #EQ 1.35, lagged circulatory normal force 

    # @show Cpotn, Dp #both of these are starting rather large. 

    # @show Cpn #It blows up, starting at the first time step, but even more so on the second. 
    states[3] = alphaf = Cpn/Cc_na + alpha0 #EQ 1.34, delayed effective angle of incidence  #Todo: This is getting fairly large on the third time step. #Todo: I think that this C^c_{n\alpha}(s,M) should be a circulatory normal coefficient. -> If it is the compressibility corrected slope... that makes sense. But ... bummer I had a question but now it's gone. 

    states[20] = fp = seperationpoint(alphaf, alpha0, alpha1, alpha2, S1, S2, S3, S4) #EQ 1.33 
    states[17] = Df = Df_m*exp(-deltas/Tf) + (fp - fp_m)*exp(-deltas/(2*Tf)) #EQ 1.36b
    states[22] = fpp = fp - Df #EQ 1.36, delayed effective seperation point. 
    


    ### 
    fterm = (1 + sqrt(fpp))/2
    # Cfsn = Cnc_naq + Cc_naq*(fterm)^2 #EQ 1.38, Normal force coefficient after accounting for separated flow from TE
    # Cfsc = Cpotc*eta*sqrt(fpp) #EQ 1.40 #(sqrt(fpp)-0.2) #Gonzalez modifications

    # @show Cnc_naq, Cc_naq, fterm #Both of the normal force coefficients are quite large. 



    states[26] = Cv = Cc_na*alphae*(1 - fterm)^2 #EQ 1.49, Normal force coefficient due to accumulated vorticity

    states[25] = Cvn = Cvn_m*exp(-deltas/Tv) + (Cv - Cv_m)*exp(-deltas/(2*Tv)) #EQ 1.47
    #Note: Cv has to be the same sign as Cfs_n

    # xbarcp = 0.2 #Todo: Decide if this goes in the airfoil or in the model. Also find out what it is. 
    # xvcp = xbarcp*(1-cos(pi*tauv/(Tvl)))
    # Cvm = -xvcp*Cvn

    ######## Update "other states"
    ### Test for Leading edge separation
    if Cpn > Cn1 
        LESF = true #LE separation can occur
    else
        LESF = false #Reattachment can occur
    end

    ### Test for Trailing edge separation
    if fpp < fpp_m
        TESF = true #TE separation in progress
    else
        TESF = false
    end

    ### Test for vortex advection
    if 0< tauv <= 2*Tvl
        VRTX = true #Vortex advection in progress
    else
        VRTX = false #Vortex is in wake. 
    end

    ### Vortex position reset
    if (tauv >= 1 + Tsh/Tvl) & (LESF = true)
        tauv = 0
    end

    if LESF
        tauv += deltat*2*U/c #No way given to update tauv. Doing a simple Euler step. 
    end

    ### Update sigma 1 - Tf modifications
    Delta_alpha0 = alpha - alpha0
    if TESF == true #(Separation)
        if Ka*Delta_alpha0 < 0
            sigma_1 = 2 #(Accelerate separation point movement)
        else
            if LESF == false
                sigma_1 = 1 #(LE separation can occur)
            else
                if fpp_m <= 0.7 
                    sigma_1 = 2 #(accelerate separation point movement if separation is occuring )
                else
                    sigma_1 = 1.75
                end
            end
        end
    else #(reattachment (`TESF==false`))
        if LESF == false
            sigma_1 = 0.5 #(Slow down reattachment)
        elseif VRTX == true  && 0 <= tauv <= Tvl
            sigma_1 = 0.25 #(No flow reattachment if vortex shedding is in progress)
        elseif Ka*Delta_alpha0 > 0
            sigma_1 = 0.75
        end 
    end

    ### Update sigma3 - Tv modifications
    if Tvl <= tauv <= 2*Tvl
        sigma3 = 3 #Postshedding
        if TESF== false
            sigma3 = 4 #Accelerate vortex lift decay
            if VRTX == true && 0 <= tauv <= Tvl
                if Ka*Delta_alpha0 < 0
                    sigma3 = 2 #Accelerate vortex lift decay
                else
                    sigma3 = 1 #default
                end
            end
        end
    else
        if Ka*Delta_alpha0 < 0
            sigma3 = 4 #vortex lift must decay fast
        end
    end
    
    if TESF == false && Kq*Delta_alpha0 < 0 #Note that it's Kq and not K_alpha
        sigma3 = 1 #Default
    end

    states[24] = tauv
    states[28] = sigma1
    states[29] = sigma3


    return states
end

function BLAD_coefficients(dsmodel::BeddoesLeishman, states, U, c, Cnfit, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tvl, xcp, a)

    _, A5, b5, _, eta = dsmodel.constants 

    Ka = states[6]
    Kpa = states[12]
    Kq = states[7]
    Kpq = states[13]

    M = U/a # Mach number
    beta = sqrt(1 - M^2) #Prandtl-Glauert compressibility correction factor
    TI = c/a




    ##### Prepare inputs for normal force coefficient. 
    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2)  
    k_alpha = 1/bot #Equation 1.11a
    Talpha = 3*k_alpha*TI/4
    Cnc_nalpha = 4*Talpha*(Ka-Kpa)/M

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) 
    k_q = 1/bot #Equation 1.11b
    Tq = 3*k_q*TI/4
    Cnc_nq = Tq*(Kq - Kpq)/M

    Cnc_naq = Cnc_nalpha + Cnc_nq

    alpha = states[2]
    X1 = states[8]
    X2 = states[9]
    alphae = (alpha - alpha0) - X1 - X2
    Cc_na = dcndalpha/beta #EQ 1.12, Circulatory component of the normal force coefficient response to step change in alpha. 
    Cc_naq = Cc_na*alphae

    fpp = states[22]
    fterm = (1 + sqrt(fpp))/2

    Cfsn = Cnc_naq + Cc_naq*(fterm)^2 #EQ 1.38, Normal force coefficient after accounting for separated flow from TE

    Cvn = states[25]

    ### Total normal force 
    Cn = Cfsn + Cvn #EQ 1.53





    ######### Prepare inputs for chordwise force coefficient. 
    
    Cpotc = Cc_naq*tan(alphae + alpha0) #Equation 1.21  
    Cfsc = Cpotc*eta*sqrt(fpp) #EQ 1.40 #(sqrt(fpp)-0.2) #Gonzalez modifications

    ### Chordwise force 
    Cc = Cfsc + Cvn*tan(alphae) #EQ 1.55






    ### Lift and Drag
    Cl = Cn*cos(alpha) + Cc*sin(alpha)
    Cd = Cn*sin(alpha) - Cc*cos(alpha) + Cd0




    
    Cc_naq = Cc_na*alphae

    q = states[5]
    Kpppq = states[15]
    Cc_mq = -dcndalpha*(q-Kpppq)*c/(16*beta*U)

    ### Noncirculatory component of moment due to change in alpha
    Cnc_ma = -Cnfit(alpha)/4 #EQ 1.27 Note: This equation was missed in the documented algorithm.

    bot = 15*(1-M) + 3*dcndalpha*A5*b5*beta*M*M/2  
    kmq = 7/bot #EQ 1.29b
    Kppq = states[14]
    Cnc_mq = -7*TI*kmq*kmq*(Kq-Kppq)/(12*M) #EQ 1.29

    tauv = states[24]
    xbarcp = 0.2 #Todo: Is this different from xcp? 
    xvcp = xbarcp*(1-cos(pi*tauv/(Tvl)))
    Cvm = -xvcp*Cvn

    ### Moment
    Cm = Cm0 - Cc_naq*(xcp - 0.25) + Cc_mq + Cnc_ma + Cnc_mq + Cvm 


    return Cn, Cc, Cl, Cd, Cm
end

function initialize_ADO(aoavec, tvec, airfoil::Airfoil, c, a) 
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
    X1 = X2 = X3 = X4 = 0.0
    Kpa = Kpq = Kppq = Kpppq = 0.0
    Dp = Df = Dfc = 0.0
    Cpotn = airfoil.cl(aoavec[1])
    fp = fpc = fpp = fppc = 1.0
    tauv = 0.0
    Cvn = 0.0
    Cv = 0.0
    Daf = 0.0
    sigma1 = sigma3 = 1.0

    states = [aoa, alpha, alphaf, aoadot, q, Ka, Kq, X1, X2, X3, X4, Kpa, Kpq, Kppq, Kpppq, Dp, Df, Dfc, Cpotn, fp, fpc, fpp, fppc, tauv, Cvn, Cv, Daf, sigma1, sigma3]

    loads = [airfoil.cl(aoa), airfoil.cd(aoa), airfoil.cl(aoa), airfoil.cd(aoa), airfoil.cm(aoa)]


    flags = [false, false, false]

    p = [c, a, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, alpha1, alpha2, S1, S2, S3, S4, xcp]

    return states, loads, vcat(p, flags)
end

