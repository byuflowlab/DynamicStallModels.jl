#=
AeroDyn's implementation of the Beddoes-Leishman model with Gonzalez's modifications, similar to the method implemented in OpenFAST v3.3.0.

The states that the theory doc says are required can be found in the BeddoesLeishmanAeroDyn.jl file. 

The states I use are: 
    1 - aoa - unfiltered angle of attack.
    2 - alpha - Filtered angle of attack. 
    3 - q - Unfiltered pitching rate
    4 - qf - Filtered pitching rate. 
    5 - Ka
    6 - Kq
    7 - Kpa
    8 - Kpq
    9 - X1
    10 - X2
    11 - X3
    12 - X4
    13 - Npot
    14 - Kppq
    15 - Kpppq
    16 - Dp
    17 - fp
    18 - fpc
    19 - fpm
    20 - Df
    21 - Dfc
    22 - Dfm
    23 - fpp
    24 - fppc
    25 - fppm
    26 - Cv
    27 - Nv
    28 - tauv
    29 - LESF::Bool
    30 - TESF::Bool
    31 - VRTX::Bool
    32 - firstpass::Bool
=#


#AeroDyn original implementation. 
function update_states_ADG!(airfoil::Airfoil, oldstates, states, y, deltat)  #TODO: Ryan seems to think that all of this airfoil information should pass in from the airfoil struct just fine and not affect how derivatives are calculated. 

    ### Unpack airfoil constants
    dcndalpha = airfoil.dcndalpha
    alpha0 = airfoil.alpha0
    # _, alpha1 = airfoil.alphasep
    c = airfoil.c
    # xcp = airfoil.xcp

    ### Unpack model constants
    model = airfoil.model
    a = model.a
    A1, A2, A5 = model.A
    b1, b2, b5 = model.b
    Tp, Tf0, Tv0, Tvl, Tsh = model.T
    zeta = model.zeta
    Cd0 = model.Cd0
    Cm0 = model.Cm0
    Cn1 = model.Cn1

    ### Unpack environmental inputs
    U, _, aoa, _ = y

    ### Unpack states
    _, alpha_m, q_m, qf_m, Ka_m, Kq_m, Kpa_m, Kpq_m, X1_m, X2_m, X3_m, X4_m, Npot_m, Kppq_m, Kpppq_m, Dp_m, fp_m, fpc_m, fpm_m, Df_m, Dfc_m, Dfm_m, fpp_m, fppc_m, fppm_m, Cv_m, Nv_m, tauv, LESF_m, TESF_m, VRTX_m, firstpass_m = oldstates #The underscore m means that it is the previous time step (m comes before n).

    ########### Algorithm ############### (Converted from UA documentation)
    ### Initial constants
    TI = c/a # Equation 1.11
    deltas = 2*U*deltat/c# Equation 1.5b #Checked against OpenFAST v3.3.0

    states[1] = aoa #Unfiltered angle of attack. 


    ###Low-pass-filtering #TODO: Might put this in a function. lowpass()
    lp_cutoff = max(1, U)*zeta/(pi*c) #Copied from OpenFAST v3.3.0. Original equation is 1.8g. 
    Clp = exp(-2*pi*deltat*lp_cutoff)
    plC = 1 - Clp


    states[2] = alpha = Clp*alpha_m + plC*aoa #low-pass-filtered angle of attack. #Equation 1.8a

    Ka = (alpha - alpha_m)/deltat #EQ 1.7b, an intermediate calculation. The real K_alpha value is calculated below. 


    states[3] = q = Ka*c/U #EQ 1.8d
    if firstpass_m==1
        q_m = q
    end
    states[4] = qf = Clp*qf_m + plC*q #low-pass-filter pitch rate. #Equation 1.8c

    

    states[5] = Ka = qf*U/c #Equation 1.8d  #This forumulation doesn't break Patankar's rules for the time step. 

    Kq = (q - q_m)/deltat #Equation 1.8e 
    if firstpass_m==1
        Kq_m = 0
    end
    states[6] = Kq = Clp*Kq_m + plC*Kq #Equation 1.8f 






    ### Update sigma 1 - Tf modifications
    Delta_alpha0 = alpha - alpha0
    sigma1 = 1
    sigma3 = 1
    sigma1c = 1
    sigma1m = 1

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
    Tfm = Tf0/sigma1m



    ### Non-circulatory components of normal force
    states[7] = Kpa = Kpa_m*exp(-deltat/Talpha) + (Ka - Ka_m)*exp(-deltat/(2*Talpha)) #Equation 1.18b, Deficiency function for eq. 1.18b
    states[8] = Kpq = Kpq_m*exp(-deltat/Tq) + (Kq - Kq_m)*exp(-deltat/(2*Tq)) #Equation 1.19b, deficiency function for eq. 1.19b
    
    Nnoncirc_a = 4*Talpha*(Ka-Kpa)/M #Equation 1.18a, Noncirculatory component of normal force due to changes in alpha
    Nnoncirc_q = -Tq*(Kq - Kpq)/M #Equation 1.19a, Noncirculatory component of normal force due to changes in pitching rate. 

    Nnoncirc_aq = Nnoncirc_a + Nnoncirc_q #Equation 1.17, Noncirculatory component of normal force via superposition. 





    ### Update states 1 and 2
    beta2 = beta^2 #beta squared
    delta_alpha = alpha - alpha_m #Note: Not explicitly defined in the documentation. Assumed. 
    deltaq = qf - qf_m #Not explicityly stated in the docs. Assumed

    states[9] = X1 = X1_m*exp(-b1*beta2*deltas) + A1*exp(-b1*beta2*deltas/2)*delta_alpha #EQ 1.15a
    states[10] = X2 = X2_m*exp(-b2*beta2*deltas) + A2*exp(-b2*beta2*deltas/2)*delta_alpha #EQ 1.15b
    states[11] = X3 = X3_m*exp(-b1*beta2*deltas) + A1*exp(-b1*beta2*deltas/2)*deltaq #EQ 1.16a #Checked against OpenFAST v3.3.0 #Note: All other models set X3 and X4 and Cn_q_circ to zero.
    states[12] = X4 = X4_m*exp(-b2*beta2*deltas) + A2*exp(-b2*beta2*deltas/2)*deltaq #EQ 1.16a #Diddo

    ### Circulatory component of normal force
    alphae = (alpha - alpha0) - X1 - X2 #EQ 1.14, Effective angle of attack
    dcndalpha_circ = dcndalpha/beta #EQ 1.12, Circulatory component of the normal force coefficient response to step change in alpha. #Checked against OpenFAST v3.3.0

    # @show dcndalpha, beta

    Ncirc_aq = dcndalpha_circ*alphae #+ Cc_nq #EQ 1.13, Circulatory component of normal force via lumped approach. #TODO: This appears to be equal to Cpotcn

    ### Total normal force under attached conditions
    states[13] = Npot = Ncirc_aq + Nnoncirc_aq #Equation 1.20




    ### Noncirculatory component of moment due to change in pitch 
    bot = 15*(1-M) + 3*dcndalpha*A5*b5*beta*M*M/2  
    kmq = 7/bot #EQ 1.29b #Checked that uses slope. 
    states[14] = Kppq = Kppq_m*exp(-deltat/(kmq*kmq*TI)) + (Kq - Kq_m)*exp(-deltat/(2*kmq*kmq*TI)) #EQ 1.29c

    ### Circulatory component of moment.  
    states[15] = Kpppq = Kpppq_m*exp(-b5*beta2*deltas) + A5*deltaq*exp(-b5*beta2*deltas/2) #EQ 1.26 





    ####### boundary layer response 
    if firstpass_m==1
        Npot_m = Npot
    end
    states[16] = Dp = Dp_m*exp(-deltas/Tp) + (Npot - Npot_m)*exp(-deltas/(Tp*2)) #EQ1.35b, deficiency function. 

    Cpn = Npot - Dp #EQ 1.35, lagged circulatory normal force 

    # @show Npot, Dp


    alphaf = Cpn/dcndalpha_circ + alpha0 #EQ 1.34, delayed effective angle of incidence  Note: this is not the filtered angle of attack. 

    # @show alphaf, Cpn, dcndalpha_circ, alpha0

    
    #### Separation points
    states[17] = fp = separationpoint(airfoil, alphaf, dcndalpha_circ) #EQ 1.33, modifications from OpenFAST v3.3.0
    states[18] = fpc = chordwiseseparationpoint(airfoil, alphaf, dcndalpha_circ) 

    Cntemp = airfoil.cl(alphaf)*cos(alphaf) + (airfoil.cd(alphaf) - Cd0)*sin(alphaf) #TODO: Do I want too move this into a function? Like the who thing to get the moment separation point value. 

    if abs(Cntemp)<0.01
        states[19] = fpm = 0 #TODO: Theoretically this makes more sense to be a one, but... whatever this is what OpenFAST v3.3.0 has. 
    else
        states[19] = fpm = (airfoil.cm(alphaf)-Cm0)/Cntemp
    end

    


    #### Deficiency functions
    if firstpass_m==1
        states[20] = Df = 0.0
    else
        states[20] = Df = Df_m*exp(-deltas/Tf) + (fp - fp_m)*exp(-deltas/(2*Tf)) #EQ 1.36b
    end

    
    if firstpass_m==1
        states[21] = Dfc = 0.0
    else
        states[21] = Dfc = Dfc_m*exp(-deltas/Tfc) + (fpc - fpc_m)*exp(-deltas/(2*Tfc)) #EQ 1.36b Chord wise separation point deficiency function
    end

    if firstpass_m == 1.0 
        states[22] = Dfm = 0
    else
    
        states[22] = Dfm = Dfm_m*exp(-deltas/Tfm) + (fpm - fpm_m)*exp(-deltas/(2*Tfm)) # Moment separation point deficiency function
    end




    #### Delayed Effective separation points
    states[23] = fpp = fp - Df #EQ 1.36, normal delayed effective seperation point. 
    states[24] = fppc = fpc - Dfc #EQ 1.36a chordwise
    states[25] = fppm = fpm - Dfm # Moment. 



    fterm = (1 + 2*sqrt(fpp))/3
    states[26] = Cv = dcndalpha_circ*alphae*(1 - fterm)^2 #EQ 1.50, Normal force coefficient due to accumulated vorticity 

    
    
    if firstpass_m==1
        states[27] = Nv = 0.0
    else
        if (tauv>Tvl)&&(Ka*Delta_alpha0>0) #TODO: I wonder if this get's taken care of by the "other states" updates. -> I don't know what I meant by that anymore. 
            states[27] = Nv = Nv_m*exp(-2*deltas/Tv)
        else
            states[27] = Nv = Nv_m*exp(-deltas/Tv) + (Cv - Cv_m)*exp(-deltas/(2*Tv)) #EQ 1.47
        end
    end

    if Nv<0 #TODO: Is this going to be a problem for optimization? 
        states[27] = 0
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

    

    states[28] = tauv
    states[29] = LESF
    states[30] = TESF
    states[31] = VRTX
    states[32] = firstpass = 0.0
end



function BLADG_coefficients(airfoil::Airfoil, states, y) #Todo: I don't know if I need this function. 
    loads = zeros(3)
    BLADG_coefficients!(airfoil, loads, states, y)
    return loads
end


function BLADG_coefficients!(airfoil::Airfoil, loads, states, y)

    U = y[1]

    c = airfoil.c
    clfit = airfoil.cl
    cdfit = airfoil.cd
    dcndalpha = airfoil.dcndalpha
    alpha0 = airfoil.alpha0
    xcp = airfoil.xcp

    model = airfoil.model

    Cd0 = model.Cd0
    Cm0 = model.Cm0
    A1, A2, A5 = model.A
    b1, b2, b5 = model.b
    Tvl = model.T[4]
    eta = model.eta
    a = model.a

    aoa = states[1]
    alpha = states[2]

    qf = states[4]
    Ka = states[5]
    Kq = states[6]
    Kpa = states[7]
    Kpq = states[8]
    X1 = states[9]
    X2 = states[10]
    X3 = states[11] 
    X4 = states[12]

    Kppq = states[14]
    Kpppq = states[15]

    fpp = states[23]
    fppc = states[24]
    fppm = states[25]

    Cvn = states[27]
    tauv = states[28]

    firstpass = states[32]


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




    alphae = (alpha - alpha0) - X1 - X2 #EQ 1.14, Effective angle of attack

    dcndalpha_circ = dcndalpha/beta #EQ 1.12, Circulatory component of the normal force coefficient response to step change in alpha. #Checked that uses slope. 
    # Ncirc_q = dcndalpha_circ*(q - X3 - X4)/2 #EQ 1.16 
    Ncirc_q = dcndalpha_circ*qf/2 - X3 -X4 #EQ 1.16 as Openfast v3.3.0 shows it. #TODO: This can't be right. It's asking for the circulatory normal force due to pitching, and it's using states X3 and X4.. I thought those states were for the chordwise force.... 

    Ncirc_aq = dcndalpha_circ*alphae #+ Cc_nq #EQ 1.13, Circulatory component of normal force via lumped approach.

    

    
    fterm = (1 + 2*sqrt(fpp))/3

    Cfsn = Nnoncirc_aq + Ncirc_q + Ncirc_aq*(fterm^2) #Equation from OpenFAST v3.3.0. From Equations 1.38 & 1.39



    ### Total normal force 
    if tauv>0 #OpenFAST v3.3.0 - UnsteadyAero.f90 line 3230  
        Cn = Cfsn + Cvn #EQ 1.53b
    else
        Cn = Cfsn
    end



    ######### Prepare inputs for chordwise force coefficient. 
    # Cpot = Ncirc_aq*tan(alphae + alpha0) #Equation 1.21  
    Cpot = dcndalpha_circ*alphae*aoa #OpenFAST v3.3.0 #They change to this modification because they say it performs better. 

    
    Cfsc = Cpot*eta*(sqrt(fppc)-0.2) #EQ 1.40, Gonzalez modifications 

    if c==2.086
        # @show Cpot, fppc
    end

    ### Chordwise force 
    Cc = Cfsc #EQ 1.55b




    ### Lift and Drag
    Cl = Cn*cos(aoa) + Cc*sin(aoa)
    Cd = Cn*sin(aoa) - Cc*cos(aoa) + Cd0 #Adding frictional drag back in. 

    
    if firstpass==1.0
        Cl = clfit(aoa)
        Cd = cdfit(aoa)
    end

    # ### Add viscoousity back into the normal and tangent coefficients
    # Cn = Cl*cos(aoa) + Cd*sin(aoa) #Cd has friction added back in. 
    # Cc = Cl*sin(aoa) - Cd*cos(aoa)


    ### Moment
    Mcirc_q = -dcndalpha*(qf-Kpppq)*c/(16*beta*U) #Equation 1.22c

    
    Mnoncirc_alpha = -Nnoncirc_a/4 #EQ 1.27 Noncirculatory component of the normal force coefficient response to step change in alpha

    bot = 15*(1-M) + 3*dcndalpha*A5*b5*beta*M*M/2  
    kmq = 7/bot #EQ 1.29b #Checked that uses slope. 
    Mnoncirc_q = -7*TI*kmq*kmq*(Kq-Kppq)/(12*M) #EQ 1.29, Circulatory component of moment.


    Cm_common = Mcirc_q + Mnoncirc_alpha + Mnoncirc_q
    Mfs = Cm0 + Cfsn*fppm + Cm_common

    Mv = xcp*(1-cos(pi*tauv/Tvl))*Cvn #EQ 1.57


    if tauv<=0
        Cm = Mfs
    else
        Cm = Mfs + Mv
    end

    loads[:] = [Cl, Cd, Cm]
    # return Cn, Cc, Cl, Cd, Cm
end

function initialize_ADG(airfoil::Airfoil, tvec, y) 

    model = airfoil.model
    Cd0 = model.Cd0

    alpha = y[3] #No delay to begin. 
    # alphadot = y[4]

    # Cn1 = airfoil.cl(alpha1)*cos(alpha1) + (airfoil.cd(alpha1) - Cd0)*sin(alpha1) #Removing the effect of frictional drag

    aoa = alpha
    # q = qf = alphadot
    aoadot = q = qf = Ka = 0.0

    Kq = 0.0
    X1 = X2 = X3 = X4 = 0.0
    Kpa = Kpq = Kppq = Kpppq = 0.0
    Dp = Df = Dfc = Dfm = 0.0
    
    Cpotn = airfoil.cl(alpha) +(airfoil.cd(alpha)-Cd0)*sin(alpha)
    fp = fpp = 1.0
    # fp = fpp = 0.0
    fpc = fppc = fclimit
    # fpc = fppc = 0.0
    fpm = fppm = 0.0 #TODO: I'm not sure what to start this as. -> It works as is, so I'll leave it. 
    tauv = 0.0
    Nv = 0.0
    Cv = 0.0

    LESF = TESF = VRTX = 0
    firstpass = 1

    states = [aoa, alpha, q, qf, Ka, Kq,  Kpa, Kpq, X1, X2, X3, X4, Cpotn, Kppq, Kpppq, Dp, fp, fpc, fpm, Df, Dfc, Dfm, fpp, fppc, fppm, Cv, Nv, tauv, LESF, TESF, VRTX, firstpass]


    loads = [airfoil.cl(aoa), airfoil.cd(aoa), airfoil.cm(aoa)]


    return states, loads
end




##########################################################################################
################################ Functions for Testing ###################################
##########################################################################################


function getintermediatestates(states, U, a, c, dcndalpha, A1, A2, b1, b2, alpha0)
    Ka = states[5]
    Kq = states[6]
    Kpa = states[7]
    Kpq = states[8]

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

    alpha = states[2]
    qf = states[4]

    X1 = states[9]
    X2 = states[10]
    X3 = states[11] #EQ 1.16a 
    X4 = states[12] #EQ 1.16a
    alphae = (alpha - alpha0) - X1 - X2

    dcndalpha_circ = dcndalpha/beta
    Ncirc_q = dcndalpha_circ*qf/2 - X3 -X4
    Ncirc_aq = dcndalpha_circ*alphae 

    fpp = states[23]
    fterm = (1 + 2*sqrt(fpp))/3

    Cfsn = Nnoncirc_aq + Ncirc_q + Ncirc_aq*(fterm^2) 

    Cvn = states[27]

    aoa = states[1]
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
        model = airfoil.model
        Cfsn[i], Cvn[i], Nnc_aq[i], Nc_q[i], Nc_aq[i], fpp[i], Nnc_a[i], Nnc_q[i], Talpha[i], M[i], k_alpha[i], TI[i], alphae[i], k_q[i], Cpot[i], dcndalpha_circ[i] = getintermediatestates(states[i,:], Uvec[i], a, c, airfoil.dcndalpha, model.A[1], model.A[2], model.b[1], model.b[2], airfoil.alpha0)
    end

    return Cfsn, Cvn, Nnc_aq, Nc_q, Nc_aq, Nnc_a, Nnc_q, Talpha, M, k_alpha, TI, alphae, k_q, Cpot, dcndalpha_circ
end