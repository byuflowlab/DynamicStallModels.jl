#=
AeroDyn's implementation of the Beddoes-Leishman model with Gonzalez's modifications, similar to the method implemented in OpenFAST v3.3.0.

=#


function getloads_BLAG(dsmodel::BeddoesLeishman, states, p, airfoil)
    c, a, dcndalpha, alpha0, Cd0, Cm0, A1, A2, A5, b1, b2, b5, Tf0, Tv0, Tp, Tvl, Tsh, Cn1, xcp, zeta, U, _ = p
    clfit = airfoil.cl #TODO: We might as well do this inside the blag coefficients. 
    cdfit = airfoil.cd
    eta = airfoil.eta

    Cn, Cc, Cl, Cd, Cm = BLADG_coefficients(dsmodel::BeddoesLeishman, states, U, c, clfit, cdfit, dcndalpha, alpha0, Cd0, Cm0, A1, A2, A5, b1, b2, b5, Tvl, xcp, eta, a)
    return [Cn, Cc, Cl, Cd, Cm]
end




#AeroDyn original implementation. 
function update_states_ADG(dsmodel::BeddoesLeishman, oldstates, c, a, U, deltat, aoa, dcndalpha, alpha0, A1, A2, A5, b1, b2, b5, Tf0, Tv0, Tp, Tvl, Tsh, Cn1, zeta, afidx)  #TODO: Ryan seems tot think that all of this airfoil information should pass in from the airfoil struct just fine and not affect how derivatives are calculated. 

    ### Unpack
    airfoil = dsmodel.airfoils[afidx]

    alpha_m, alphaf_m, q_m, Ka_m, Kq_m, X1_m, X2_m, X3_m, X4_m, Kpa_m, Kpq_m, Kppq_m, Kpppq_m, Dp_m, Df_m, Dfc_m, Npot_m, fp_m, fpp_m, fpc_m, fppc_m, tauv, Nv_m, Cv_m, LESF_m, TESF_m, VRTX_m, firstpass_m, qf_m = oldstates #The underscore m means that it is the previous time step (m comes before n).

    


    #= States
    1 - alpha # Filtered angle of attack. 
    2 - alphaf #TODO: I'm not sure that this state is needed. Note: this is not the filtered angle of attack. 
    3 - q # Unfiltered pitching rate
    4 - Ka
    5 - Kq
    6 - X1
    7 - X2
    8 - X3
    9 - X4
    10 - Kpa
    11 - Kpq
    12 - Kppq
    13 - Kpppq
    14 - Dp
    15 - Df
    16 - Dfc
    17 - Cpotn
    18 - fp
    19 - fpp
    20 - fpc
    21 - fppc
    22 - tauv
    23 - Cvn
    24 - Cv
    25 - LESF::Bool
    26 - TESF::Bool
    27 - VRTX::Bool
    28 - firstpass::Bool
    29 - qf - Filtered pitching rate. 
    30 - aoa - unfiltered angle of attack. 
    =#


    states = zeros(30) #TODO: Consider putting a function to return this value. 

    states[30] = aoa #Unfiltered angle of attack. 


    

    ########### Algorithm ############### (Converted from UA documentation)
    ### Initial constants
    TI = c/a # Equation 1.11
    deltas = 2*U*deltat/c# Equation 1.5b #Checked against OpenFAST v3.3.0

    #TODO: I don't think that I need to set alpha_m because I've already set it. 
    if firstpass_m==1.0
        alpha_m = aoa
        alpha_f_m = aoa
    end


    ###Low-pass-filtering #TODO: Might put this in a function. lowpass()
    lp_cutoff = max(1, U)*zeta/(pi*c) #Copied from OpenFAST v3.3.0. Original equation is 1.8g. 
    Clp = exp(-2*pi*deltat*lp_cutoff)
    plC = 1 - Clp


    states[1] = alpha = Clp*alpha_m + plC*aoa #low-pass-filtered angle of attack. #Equation 1.8a

    Ka = (alpha - alpha_m)/deltat #EQ 1.7b
    if firstpass_m==1 #TODO. Does this even get used? I'm not storing Ka... so I don't think this gets used... which means I probably have a misnaming on Ka and Kalpha... great -> But I'm getting q right, so I don't know if it matters. It gets used below when I update the "other states". 
        Ka_m = 0
    end




    states[3] = q = Ka*c/U #EQ 1.8d
    if firstpass_m==1
        q_m = q
    end
    states[29] = qf = Clp*qf_m + plC*q #low-pass-filter pitch rate. #Equation 1.8c

    

    states[4] = Ka = qf*U/c #Equation 1.8d  #This forumulation doesn't break Patankar's rules for the time step. 

    Kq = (q - q_m)/deltat #Equation 1.8e 
    if firstpass_m==1
        Kq_m = 0
    end
    states[5] = Kq = Clp*Kq_m + plC*Kq #Equation 1.8f 





    ### Update sigma 1 - Tf modifications
    Delta_alpha0 = alpha - alpha0
    sigma1 = 1
    sigma3 = 1
    sigma1c = 1

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

    ### Update sigma 1c - Tfc modifications -> It appears that in OpenFASTv3.3.0 - UnsteadyAero.f90 they never update sigma1c from 1.0... so... it remains constant? 
    # if TESF_m == 1 #(Separation)
    #     if Ka_m*Delta_alpha0 < 0
    #         sigma1c = 2 #(Accelerate separation point movement)
    #     else
    #         if LESF_m == 0
    #             sigma1c = 1 #(LE separation can occur)
    #         else
    #             if fppc_m <= 0.7 
    #                 sigma1c = 2 #(accelerate separation point movement if separation is occuring )
    #             else
    #                 sigma1c = 1.75
    #             end
    #         end
    #     end
    # else #(reattachment (`TESF==false`))
    #     if LESF_m == 0
    #         sigma1c = 0.5 #(Slow down reattachment)
    #     elseif VRTX_m == 1  && 0 <= tauv <= Tvl
    #         sigma1c = 0.25 #(No flow reattachment if vortex shedding is in progress)
    #     elseif Ka_m*Delta_alpha0 > 0
    #         sigma1c = 0.75
    #     end 
    # end



    ### More constants
    M = U/a # Mach number
    beta = sqrt(1 - M^2) #Prandtl-Glauert compressibility correction factor #TODO: Need something to cap if M>1. Or at least some sort of diversion behavior. 

    bot = (1-M) + dcndalpha*(M^2)*beta*(A1*b1 + A2*b2)/2 #TODO: Calculated solely at first time step.... Maybe we pull this out and pass it in as a function argument. ... Or just calculate it every time as it isn't super costly. 
    k_alpha = 1/bot #Equation 1.11a #Checked that uses slope. 

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) #TODO. The documentation has the same equation for 1.11b as 1.11a. -> In Leishman's 1900 state space paper, he has the second term in the denominator of k_alpha to be half of what is shown here. I'll try what I have, and if it is off, then I'll try changing it. -> It appears to be correct. -> OpenFAST v3.3.0 has a one half on kalpha. It made very little difference.  #Todo/ Double check that both of these are correct. I'm pretty sure they are, but double check cause my note conflicts with what I have. -> The 1/2 is on the denominator of k_alpha in OpenFAST. 
    k_q = 1/bot #Equation 1.11b #Checked that uses slope. 

    Talpha = 3*k_alpha*TI/4 #Equation 1.10a
    Tq = 3*k_q*TI/4 #Equation 1.10b

    Tf = Tf0/sigma1 #Equation 1.37
    Tv = Tv0/sigma3 #Equation 1.48

    Tfc = Tf0/sigma1c



    ### Non-circulatory components of normal force
    states[10] = Kpa = Kpa_m*exp(-deltat/Talpha) + (Ka - Ka_m)*exp(-deltat/(2*Talpha)) #Equation 1.18b, Deficiency function for eq. 1.18b
    Nnoncirc_a = 4*Talpha*(Ka-Kpa)/M #Equation 1.18a, Noncirculatory component of normal force due to changes in alpha

    states[11] = Kpq = Kpq_m*exp(-deltat/Tq) + (Kq - Kq_m)*exp(-deltat/(2*Tq)) #Equation 1.19b, deficiency function for eq. 1.19b
    
    Nnoncirc_q = -Tq*(Kq - Kpq)/M #Equation 1.19a, Noncirculatory component of normal force due to changes in pitching rate. 

    Nnoncirc_aq = Nnoncirc_a + Nnoncirc_q #Equation 1.17, Noncirculatory component of normal force via superposition. 

    




    ### Update states 1 and 2
    beta2 = beta^2 #beta squared
    delta_alpha = alpha - alpha_m #Note: Not explicitly defined in the documentation. Assumed. 
    deltaq = qf - qf_m #Not explicityly stated in the docs. Assumed

    states[6] = X1 = X1_m*exp(-b1*beta2*deltas) + A1*exp(-b1*beta2*deltas/2)*delta_alpha #EQ 1.15a
    states[7] = X2 = X2_m*exp(-b2*beta2*deltas) + A2*exp(-b2*beta2*deltas/2)*delta_alpha #EQ 1.15b



    ### Circulatory component of normal force
    alphae = (alpha - alpha0) - X1 - X2 #EQ 1.14, Effective angle of attack

    dcndalpha_circ = dcndalpha/beta #EQ 1.12, Circulatory component of the normal force coefficient response to step change in alpha. #Checked against OpenFAST v3.3.0

    #Note: All other models set X3 and X4 and Cn_q_circ to zero. 
    states[8] = X3 = X3_m*exp(-b1*beta2*deltas) + A1*exp(-b1*beta2*deltas/2)*deltaq #EQ 1.16a #Checked against OpenFAST v3.3.0
    states[9] = X4 = X4_m*exp(-b2*beta2*deltas) + A2*exp(-b2*beta2*deltas/2)*deltaq #EQ 1.16a #Diddo

    Ncirc_aq = dcndalpha_circ*alphae #+ Cc_nq #EQ 1.13, Circulatory component of normal force via lumped approach. #TODO: This appears to be equal to Cpotcn

    

    ### Circulatory component of moment. #Question: Why did they calculate the moment here? Do I need some of these things here? Is there a better spot to put this? 
    states[13] = Kpppq = Kpppq_m*exp(-b5*beta2*deltas) + A5*deltaq*exp(-b5*beta2*deltas/2) #EQ 1.26



    ### Total normal force under attached conditions
    states[17] = Npot = Ncirc_aq + Nnoncirc_aq #Equation 1.20
    


    ### Noncirculatory component of moment due to change in pitch #Note: If this can, this should probably go with the other moment calculation. 
    bot = 15*(1-M) + 3*dcndalpha*A5*b5*beta*M*M/2  
    kmq = 7/bot #EQ 1.29b #Checked that uses slope. 
    states[12] = Kppq = Kppq_m*exp(-deltat/(kmq*kmq*TI)) + (Kq - Kq_m)*exp(-deltat/(2*kmq*kmq*TI)) #EQ 1.29c



    ####### boundary layer response 
    if firstpass_m==1
        Npot_m=Npot
    end
    states[14] = Dp = Dp_m*exp(-deltas/Tp) + (Npot - Npot_m)*exp(-deltas/(Tp*2)) #EQ1.35b, deficiency function. 
    Cpn = Npot - Dp #EQ 1.35, lagged circulatory normal force 


    states[2] = alphaf = Cpn/dcndalpha_circ + alpha0 #EQ 1.34, delayed effective angle of incidence  

    
    states[18] = fp = separationpoint(airfoil, alphaf, dcndalpha_circ) #EQ 1.33, modifications from OpenFAST v3.3.0


    if firstpass_m==1
        states[15] = Df = 0.0
    else
        states[15] = Df = Df_m*exp(-deltas/Tf) + (fp - fp_m)*exp(-deltas/(2*Tf)) #EQ 1.36b
    end

    states[19] = fpp = fp - Df #EQ 1.36, delayed effective seperation point. 


    states[20] = fpc = chordwiseseparationpoint(airfoil, alphaf, dcndalpha_circ) 

    if firstpass_m==1
        states[16] = Dfc = 0.0
    else
        states[16] = Dfc = Dfc_m*exp(-deltas/Tfc) + (fpc - fpc_m)*exp(-deltas/(2*Tfc)) #EQ 1.36b (applied to the chordwise force)
    end

    states[21] = fppc = fpc - Dfc #EQ 1.36a (applied to the chordwise force). 
    



    fterm = (1 + 2*sqrt(fpp))/3

    states[24] = Cv = dcndalpha_circ*alphae*(1 - fterm)^2 #EQ 1.50, Normal force coefficient due to accumulated vorticity 

    # states[23] = Nv = Nv_m*exp(-deltas/Tv) + (Cv - Cv_m)*exp(-deltas/(2*Tv)) #EQ 1.47 #TODO: They have something different. See below. 
    if firstpass_m==1
        states[23] = Nv = 0.0
    else
        if (tauv>Tvl)&&(Ka*Delta_alpha0>0) #TODO: I wonder if this get's taken care of by the "other states" updates
            states[23] = Nv = Nv_m*exp(-2*deltas/Tv)
        else
            states[23] = Nv = Nv_m*exp(-deltas/Tv) + (Cv - Cv_m)*exp(-deltas/(2*Tv))
        end
    end

    if states[23]<0 #TODO: Is this going to be a problem for optimization? 
        states[23] = 0
    end
    


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

    

    states[22] = tauv
    states[25] = LESF
    states[26] = TESF
    states[27] = VRTX
    states[28] = firstpass = 0.0



    return states
end

function BLADG_coefficients(dsmodel::BeddoesLeishman, states, U, c, af::Airfoil, a)
    clfit = af.cl
    cdfit = af.cd
    dcndalpha = af.dcndalpha
    alpha0 = af.alpha0
    Cd0 = af.cd(alpha0)
    Cm0 = af.cm(alpha0)
    A1 = af.A[1]
    A2 = af.A[2]
    A5 = af.A[3]
    b1 = af.b[1]
    b2 = af.b[2]
    b5 = af.b[3]
    Tvl = af.T[4]

    xcp = af.xcp
    eta = af.eta

    return BLADG_coefficients(dsmodel, states, U, c, clfit, cdfit, dcndalpha, alpha0, Cd0, Cm0, A1, A2, A5, b1, b2, b5, Tvl, xcp, eta, a)
end

function BLADG_coefficients(dsmodel::BeddoesLeishman, states, U, c, clfit, cdfit, dcndalpha, alpha0, Cd0, Cm0, A1, A2, A5, b1, b2, b5, Tvl, xcp, eta, a)

    # _, A5, b5, _ = dsmodel.constants 

    Ka = states[4]
    Kpa = states[10]
    Kq = states[5]
    Kpq = states[11]

    M = U/a # Mach number
    beta = sqrt(1 - M^2) #Prandtl-Glauert compressibility correction factor
    TI = c/a




    ##### Prepare inputs for normal force coefficient. 
    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2)/2
    k_alpha = 1/bot #Equation 1.11a #Checked that uses slope. 
    Talpha = 3*k_alpha*TI/4 #Equation 1.10a
    Nnoncirc_a = 4*Talpha*(Ka-Kpa)/M #Equation 1.18, Noncirculatory component of normal force due to changes in alpha


    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) 
    k_q = 1/bot #Equation 1.11b #Checked that uses slope. 
    Tq = 3*k_q*TI/4 #Equation 1.10b
    Nnoncirc_q = -Tq*(Kq - Kpq)/M #Equation 1.19

    Nnoncirc_aq = Nnoncirc_a + Nnoncirc_q #Equation 1.17, Noncirculatory component of normal force via superposition.



    alpha = states[1]
    qf = states[29]

    X1 = states[6]
    X2 = states[7]
    X3 = states[8] 
    X4 = states[9]
    alphae = (alpha - alpha0) - X1 - X2 #EQ 1.14, Effective angle of attack

    dcndalpha_circ = dcndalpha/beta #EQ 1.12, Circulatory component of the normal force coefficient response to step change in alpha. #Checked that uses slope. 
    # Ncirc_q = dcndalpha_circ*(q - X3 - X4)/2 #EQ 1.16 
    Ncirc_q = dcndalpha_circ*qf/2 - X3 -X4 #EQ 1.16 as Openfast v3.3.0 shows it. #TODO: This can't be right. It's asking for the circulatory normal force due to pitching, and it's using states X3 and X4.. I thought those states were for the chordwise force.... 

    Ncirc_aq = dcndalpha_circ*alphae #+ Cc_nq #EQ 1.13, Circulatory component of normal force via lumped approach.

    

    fpp = states[19]
    
    fterm = (1 + 2*sqrt(fpp))/3

    # Cfsn = Cnc_naq + Cc_naq*(fterm)^2 #+ Cc_nq #EQ 1.38, Normal force coefficient after accounting for separated flow from TE #TODO. I'm not sure that this matches Gonzalez's modifcations. #Note. There is a typo in their equation. The first part of 1.53b says that it is Cfsn, and the second part is missing the Cc_nq term from 1.39, thus either they added an extra term in 1.39, or they missed a term in 1.53b. -> It didn't make much of a difference in the NREL 5MW verification case. -> Todo. I need to check OpenFAST and see what they do. The equation I have below is the equation that OpenFAST v3.3.0 has for Cfsn. 
    # Cfsn = Nnoncirc_aq + Ncirc_q + dcndalpha_circ*alphae*(fterm^2) #EQ 1.39
    Cfsn = Nnoncirc_aq + Ncirc_q + Ncirc_aq*(fterm^2) #Equation from OpenFAST v3.3.0

    Cvn = states[23]
    tauv = states[22]

    ### Total normal force 
    if tauv>0 #OpenFAST v3.3.0 - UnsteadyAero.f90 line 3230  
        Cn = Cfsn + Cvn #EQ 1.53b
    else
        Cn = Cfsn
    end



    ######### Prepare inputs for chordwise force coefficient. 
    # Cpot = Ncirc_aq*tan(alphae + alpha0) #Equation 1.21   #Todo. What is Npot,circ? -> It's Ncirc_aq
    aoa = states[30]
    Cpot = dcndalpha_circ*alphae*aoa #OpenFAST v3.3.0 # At first I thought I had misread... but no... they definitely have C_nalpha_circ.... which suggests that maybe they mistyped. Or something? It seems like an odd correction, but it's what they have. 

    fppc = states[21]
    
    Cfsc = Cpot*eta*(sqrt(fppc)-0.2) #EQ 1.40, Gonzalez modifications 

    ### Chordwise force 
    Cc = Cfsc #EQ 1.55b




    ### Lift and Drag
    Cl = Cn*cos(aoa) + Cc*sin(aoa)
    Cd = Cn*sin(aoa) - Cc*cos(aoa) + Cd0 #Adding frictional drag back in. 

    firstpass = states[28]
    if firstpass==1.0
        Cl = clfit(aoa)
        Cd = cdfit(aoa)
    end



    # q = states[3]
    # Kpppq = states[13]
    # Cc_mq = -dcndalpha*(q-Kpppq)*c/(16*beta*U) #Equation 1.22c #Checked that uses slope. 

    ### Noncirculatory component of moment due to change in alpha
    # Cnc_ma = -Cnc_nalpha/4 #EQ 1.27 Noncirculatory component of the normal force coefficient response to step change in alpha  Note: This equation was missed in the documented algorithm. #Todo. Should I use the dynamic normal coefficient or the static noormal coefficient. Neither. It should be using equation 1.18. :) 

    # bot = 15*(1-M) + 3*dcndalpha*A5*b5*beta*M*M/2  
    # kmq = 7/bot #EQ 1.29b #Checked that uses slope. 
    # Kppq = states[12]
    # Cnc_mq = -7*TI*kmq*kmq*(Kq-Kppq)/(12*M) #EQ 1.29, Circulatory component of moment.

    
    # xbarcp = 0.2 #Todo: Is this different from xcp? They call it x bar bar cp
    # xvcp = xbarcp*(1-cos(pi*tauv/(Tvl))) #1.57b
    # Cvm = -xvcp*Cvn #1.57a

    ### Moment
    # Cm = Cn*fpp + Cc_mq + Cnc_ma + Cnc_mq + Cvm #Equation 1.60 #Todo: I don't think that this has seen the correction from Gonzalez, see equation 1.45 and 1.60

    ### Add viscoousity back into the normal and tangent coefficients
    Cn = Cl*cos(aoa) + Cd*sin(aoa) #Cd has friction added back in. 
    Cc = Cl*sin(aoa) - Cd*cos(aoa)
    Cm = 1

    return Cn, Cc, Cl, Cd, Cm
end

function initialize_ADG(Uvec, aoavec, tvec, airfoil::Airfoil, c, a) 
    # Cnfit = airfoil.cl
    dcndalpha = airfoil.dcndalpha
    alpha0 = airfoil.alpha0
    A1, A2, A5 = airfoil.A
    b1, b2, b5 = airfoil.b
    Tp, Tf0, Tv0, Tvl, Tsh = airfoil.T
    _, alpha1 = airfoil.alphasep
    xcp = airfoil.xcp
    zeta = airfoil.zeta

    Cd0 = airfoil.cd(alpha0)

    Cn1 = airfoil.cl(alpha1)*cos(alpha1) + (airfoil.cd(alpha1) - Cd0)*sin(alpha1) #Removing the effect of frictional drag
    Cm0 = airfoil.cm(alpha0)



    dt = tvec[2] - tvec[1]

    aoa = aoavec[1]
    alpha = alphaf = aoavec[1] #No delay to begin. 
    aoadot = q = Ka = 0 

    Kq = 0.0
    X1 = X2 = X3 = X4 = 0.0
    Kpa = Kpq = Kppq = Kpppq = 0.0
    Dp = Df = Dfc = 0.0
    
    Cpotn = airfoil.cl(alpha) +(airfoil.cd(alpha)-Cd0)*sin(alpha)
    fp = fpp = 1.0
    fpc = fppc = fclimit
    tauv = 0.0
    Cvn = 0.0
    Cv = 0.0

    LESF = TESF = VRTX = 0
    firstpass = 1

    states = [alpha, alphaf, q, Ka, Kq, X1, X2, X3, X4, Kpa, Kpq, Kppq, Kpppq, Dp, Df, Dfc, Cpotn, fp, fpp, fpc, fppc, tauv, Cvn, Cv, LESF, TESF, VRTX, firstpass, q, aoa]

    Cn = airfoil.cn(alpha) 
    Cc = airfoil.cc(alpha)

    loads = [Cn, Cc, airfoil.cl(aoa), airfoil.cd(aoa), airfoil.cm(aoa)]


    p = [c, a, dcndalpha, alpha0, Cd0, Cm0, A1, A2, A5, b1, b2, b5, Tf0, Tv0, Tp, Tvl, Tsh, Cn1, xcp, zeta] #20 elements (-> total of 22 elements in p. )

    envvars = [Uvec[1], aoavec[1]]

    return states, loads, vcat(p, envvars)
end

function updateenvironment_ADG(p, U, aoa) 
    p[21] = U
    p[22] = aoa
end



##########################################################################################
################################ Functions for Testing ###################################
##########################################################################################


function getintermediatestates(states, U, a, c, dcndalpha, A1, A2, b1, b2, alpha0)
    Ka = states[4]
    Kpa = states[10]
    Kq = states[5]
    Kpq = states[11]

    M = U/a # Mach number
    beta = sqrt(1 - M^2) #Prandtl-Glauert compressibility correction factor
    TI = c/a

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2)/2  
    k_alpha = 1/bot #Equation 1.11a #Checked that uses slope. 
    Talpha = 3*k_alpha*TI/4 #Equation 1.10a
    Nnoncirc_a = 4*Talpha*(Ka-Kpa)/M

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) 
    k_q = 1/bot #Equation 1.11b #Checked that uses slope. 
    Tq = 3*k_q*TI/4 #Equation 1.10b
    Nnoncirc_q = -Tq*(Kq - Kpq)/M #Equation 1.19

    Nnoncirc_aq = Nnoncirc_a + Nnoncirc_q

    alpha = states[1]
    # q = states[3]
    qf = states[29]

    X1 = states[6]
    X2 = states[7]
    X3 = states[8] #EQ 1.16a 
    X4 = states[9] #EQ 1.16a
    alphae = (alpha - alpha0) - X1 - X2

    dcndalpha_circ = dcndalpha/beta
    Ncirc_q = dcndalpha_circ*qf/2 - X3 -X4
    Ncirc_aq = dcndalpha_circ*alphae 

    fpp = states[19]
    fterm = (1 + 2*sqrt(fpp))/3

    Cfsn = Nnoncirc_aq + Ncirc_q + Ncirc_aq*(fterm^2) 

    Cvn = states[23]

    aoa = states[30]
    Cpot = dcndalpha_circ*alphae*aoa

    return Cfsn, Cvn, Nnoncirc_aq, Ncirc_q, Ncirc_aq, fpp, Nnoncirc_a, Nnoncirc_q, Talpha, M, k_alpha, TI, alphae, k_q, Cpot, dcndalpha_circ
end

function extractintermediatestates(states, Uvec, c, airfoil; a=343.0)
    n, m = size(states)

    Cfsn = zeros(n)
    Cvn = zeros(n)
    Nnc_aq = zeros(n)
    Nc_q = zeros(n)
    Nc_aq = zeros(n)
    fpp = zeros(n)
    Nnc_a = zeros(n)
    Nnc_q = zeros(n)
    Talpha = zeros(n)
    M = zeros(n)
    k_alpha = zeros(n)
    TI = zeros(n)
    alphae = zeros(n)
    k_q = zeros(n)
    Cpot = zeros(n)
    dcndalpha_circ = zeros(n)

    for i = 1:n
        Cfsn[i], Cvn[i], Nnc_aq[i], Nc_q[i], Nc_aq[i], fpp[i], Nnc_a[i], Nnc_q[i], Talpha[i], M[i], k_alpha[i], TI[i], alphae[i], k_q[i], Cpot[i], dcndalpha_circ[i] = getintermediatestates(states[i,:], Uvec[i], a, c, airfoil.dcndalpha, airfoil.A[1], airfoil.A[2], airfoil.b[1], airfoil.b[2], airfoil.alpha0)
    end

    return Cfsn, Cvn, Nnc_aq, Nc_q, Nc_aq, Nnc_a, Nnc_q, Talpha, M, k_alpha, TI, alphae, k_q, Cpot, dcndalpha_circ
end