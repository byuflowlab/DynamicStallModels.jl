#=
AeroDyn's implementation of the Beddoes-Leishman model with Gonzalez's modifications. 

=#


function getloads_BLAG(dsmodel::BeddoesLeishman, states, p, airfoil)
    c, a, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, xcp, U, _ = p
    clfit = airfoil.cl #TODO: We might as well do this inside the blag coefficients. 
    eta = airfoil.eta

    Cn, Cc, Cl, Cd, Cm = BLADG_coefficients(dsmodel::BeddoesLeishman, states, U, c, clfit, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tvl, xcp, eta, a)
    return [Cn, Cc, Cl, Cd, Cm]
end





#AeroDyn original implementation. 
function update_states_ADG(dsmodel::BeddoesLeishman, oldstates, c, a, U, deltat, aoa, dcndalpha, alpha0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, afidx)  #Todo: Ryan seems tot think that all of this airfoil information should pass in from the airfoil struct just fine and not affect how derivatives are calculated. 

    ### Unpack
    alpha_m, alphaf_m, q_m, Ka_m, Kq_m, X1_m, X2_m, X3_m, X4_m, Kpa_m, Kpq_m, Kppq_m, Kpppq_m, Dp_m, Df_m, Dfc_m, Npot_m, fp_m, fpp_m, fpc_m, fppc_m, tauv, Nv_m, Cv_m, LESF_m, TESF_m, VRTX_m = oldstates #The underscore m means that it is the previous time step (m comes before n).

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
    =#


    states = zeros(27) #Todo: Consider putting a function to return this value. 

    # @show length(states)


    zeta, A5, b5, Tsh = dsmodel.constants 
    #=
    zeta - Low-pass-fileter frequency cutoff. #TODO: Should this have the negative or should the equation have the negative? 
    Tsh - Strouhal's frequency 0.19
    =#

    ########### Algorithm ############### (Converted from UA documentation)
    ### Initial constants
    TI = c/a # Equation 1.11
    deltas = 2*U*deltat/c# Equation 1.5b #Checked against OpenFAST v3.3.0



    ###Low-pass-filtering #TODO: Might put this in a function. lowpass()
    # Clp = exp(2*pi*deltat*zeta) #low-pass-filter constant. #Equation 1.8g -> removed the negative because zeta has a negative. 
    filtercutoff = 0.5 #20 #Todo: I think this is the default. Change to pass in. -> They state in the input files that it defaults to 20, but when you print the value out it comes out as 0.5. So that's odd. 
    lp_cutoff = max(1, U)*filtercutoff/(pi*c)
    Clp = exp(-2*pi*deltat*lp_cutoff)
    plC = 1 - Clp

    # @show Clp, lp_cutoff

    states[1] = alpha = Clp*alpha_m + plC*aoa #low-pass-filtered angle of attack. #Equation 1.8a
    # @show alpha

    Ka = (alpha - alpha_m)/deltat

    q = Ka*c/U
    # q = (alpha - alpha_m)*c/(U*deltat) #Pitch rate. #Equation 1.8b #This appears to be dimensionless. 
    states[3] = q = Clp*q_m + plC*q #low-pass-filter pitch rate. #Equation 1.8c

    states[4] = Ka
    # states[4] = Ka = q*U/deltat #Equation 1.8d #Todo: this combination of q and delta t don't work for small delta t. It makes Ka super larger. So either q needs to be smaller, or delta t larger. If delta t is larger, then q should be smaller as well... which is weird, because that puts a bottom limit on delta t... which breaks CFD rules. What is it Patankar's rules? I think this breaks Patankar's rules. 
    #Note: This is Kalpha_f in OpenFAST v3.3.0
    # @show q, deltat #because deltat is small, Ka is large.... Maybe q should be smaller? 

    # if c== 1.419
    #     @show q, U, deltat
    # end

    Kq = (q - q_m)/deltat #Equation 1.8e
    states[5] = Kq = Clp*Kq_m + plC*Kq #Equation 1.8f #Todo: OpenFAST v3.3.0 uses the unfiltered q values for Kq.... which kind of is problematic since I store the filtered values. I wonder if I shouldn't apply a filter to Kq since I've applied a filter to q. Or if I should just convert to what they've done... I don't know. But I don't know if it's that big of a deal. 

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

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2)/2 #TODO: Calculated solely at first time step.... Maybe we pull this out and pass it in as a function argument. ... Or just calculate it every time as it isn't super costly. 
    k_alpha = 1/bot #Equation 1.11a #Checked that uses slope. 

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) #TODO. The documentation has the same equation for 1.11b as 1.11a. -> In Leishman's 1900 state space paper, he has the second term in the denominator of k_alpha to be half of what is shown here. I'll try what I have, and if it is off, then I'll try changing it. -> It appears to be correct. TODO: change my documentation. OpenFAST v3.3.0 has a one half on kalpha. It made very little difference. 
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
    
    Nnoncirc_q = -Tq*(Kq - Kpq)/M #Equation 1.19a, Noncirculatory component of normal force due to changes in pitching rate. #TODO. The documentation conflicts on whether or not a minus should be included here. -> I'm going to assume that it is positive. And if I'm wrong... I'll change it. It seems to be correct with positive. -> OpenFAST v3.3.0 UnsteadyAero.f90 line 496 has the negative. (No noticeable change.)

    Nnoncirc_aq = Nnoncirc_a + Nnoncirc_q #Equation 1.17, Noncirculatory component of normal force via superposition. 

    




    ### Update states 1 and 2
    beta2 = beta^2 #beta squared
    delta_alpha = alpha - alpha_m #Note: Not explicitly defined in the documentation. Assumed. 
    deltaq = q - q_m #Not explicityly stated in the docs. Assumed

    states[6] = X1 = X1_m*exp(-b1*beta2*deltas) + A1*exp(-b1*beta2*deltas/2)*delta_alpha #EQ 1.15a
    states[7] = X2 = X2_m*exp(-b2*beta2*deltas) + A2*exp(-b2*beta2*deltas/2)*delta_alpha #EQ 1.15b



    ### Circulatory component of normal force
    alphae = (alpha - alpha0) - X1 - X2 #EQ 1.14, Effective angle of attack

    dcndalpha_circ = dcndalpha/beta #EQ 1.12, Circulatory component of the normal force coefficient response to step change in alpha. #Checked against OpenFAST v3.3.0

    #Note: All other models set X3 and X4 and Cn_q_circ to zero. 
    states[8] = X3 = X3_m*exp(-b1*beta2*deltas) + A1*exp(-b1*beta2*deltas/2)*deltaq #EQ 1.16a #Checked against OpenFAST v3.3.0
    states[9] = X4 = X4_m*exp(-b2*beta2*deltas) + A2*exp(-b2*beta2*deltas/2)*deltaq #EQ 1.16a #Diddo
    Ncirc_q = dcndalpha_circ*(q - X3 - X4)/2 #EQ 1.16, Circulatory component of the normal force coefficient response to step change in pitch. #Checked against OpenFAST v3.3.0 #Todo: Where does this get used? -> In the Cn_FS calculation. 

    Ncirc_aq = dcndalpha_circ*alphae #+ Cc_nq #EQ 1.13, Circulatory component of normal force via lumped approach. #TODO: This appears to be equal to Cpotcn

    

    ### Circulatory component of moment. #Question: Why did they calculate the moment here? Do I need some of these things here? Is there a better spot to put this? 
    states[13] = Kpppq = Kpppq_m*exp(-b5*beta2*deltas) + A5*deltaq*exp(-b5*beta2*deltas/2) #EQ 1.26



    ### Total normal force under attached conditions
    states[17] = Npot = Ncirc_aq + Nnoncirc_aq #Equation 1.20
    # if c== 1.419
    #     @show Cc_naq, Cnc_naq
    # end


    ### Noncirculatory component of moment due to change in pitch #Note: If this can, this should probably go with the other moment calculation. 
    bot = 15*(1-M) + 3*dcndalpha*A5*b5*beta*M*M/2  
    kmq = 7/bot #EQ 1.29b #Checked that uses slope. 
    states[12] = Kppq = Kppq_m*exp(-deltat/(kmq*kmq*TI)) + (Kq - Kq_m)*exp(-deltat/(2*kmq*kmq*TI)) #EQ 1.29c



    ####### boundary layer response 
    states[14] = Dp = Dp_m*exp(-deltas/Tp) + (Npot - Npot_m)*exp(-deltas/(Tp*2)) #EQ1.35b, deficiency function. 
    Cpn = Npot - Dp #EQ 1.35, lagged circulatory normal force 


    states[2] = alphaf = Cpn/dcndalpha_circ + alpha0 #EQ 1.34, delayed effective angle of incidence  

    airfoil = dsmodel.airfoils[afidx]
    states[18] = fp = separationpoint(airfoil, alphaf) #EQ 1.33 
    states[15] = Df = Df_m*exp(-deltas/Tf) + (fp - fp_m)*exp(-deltas/(2*Tf)) #EQ 1.36b
    states[19] = fpp = fp - Df #EQ 1.36, delayed effective seperation point. 
    # @show fpp

    states[20] = fpc = chordwiseseparationpoint(airfoil, alphaf)
    states[16] = Dfc = Dfc_m*exp(-deltas/Tfc) + (fpc - fpc_m)*exp(-deltas/(2*Tfc)) #EQ 1.36b (applied to the chordwise force)
    states[21] = fppc = fpc - Dfc #EQ 1.36a (applied to the chordwise force). 
    
    # println("update states fpp: ", fpp, " , ", states[18])
    # @show fp, fpp, Df
    # @show fp_m, fpp_m, Df_m #Question. How is Df_m getting set to 1.0? -> Rotors was using initializeADO. 

    ### 
    fterm = (1 + 2*sqrt(fpp))/3



    states[24] = Cv = dcndalpha_circ*alphae*(1 - fterm)^2 #EQ 1.50, Normal force coefficient due to accumulated vorticity #Todo:OpenFAST has Cv = Cn_alpha_q_circ*(1-fterm)^2. I need to check that this is the same. 

    # states[23] = Nv = Nv_m*exp(-deltas/Tv) + (Cv - Cv_m)*exp(-deltas/(2*Tv)) #EQ 1.47 #Todo: They have something different. See below. 

    if (tauv>Tvl)&&(Ka*Delta_alpha0>0) #Todo: I wonder if this get's taken care of by the "other states" updates
        states[23] = Nv = Nv_m*exp(-2*deltas/Tv)
    else
        states[23] = Nv = Nv_m*exp(-deltas/Tv) + (Cv - Cv_m)*exp(-deltas/(2*Tv))
    end

    if states[23]<0
        states[23] = 0
    end

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

    

    states[22] = tauv
    states[25] = LESF
    states[26] = TESF
    states[27] = VRTX



    return states
end

function BLADG_coefficients(dsmodel::BeddoesLeishman, states, U, c, af::Airfoil, a)
    cnfit = af.cl
    dcndalpha = af.dcndalpha
    alpha0 = af.alpha0
    Cd0 = af.cd(alpha0)
    Cm0 = af.cm(alpha0)
    A1 = af.A[1]
    A2 = af.A[2]
    b1 = af.b[1]
    b2 = af.b[2]
    Tvl = af.T[4]
    xcp = af.xcp
    eta = af.eta

    return BLADG_coefficients(dsmodel, states, U, c, cnfit, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tvl, xcp, eta, a)
end

function BLADG_coefficients(dsmodel::BeddoesLeishman, states, U, c, Clfit, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tvl, xcp, eta, a)

    _, A5, b5, _ = dsmodel.constants 

    Ka = states[4]
    Kpa = states[10]
    Kq = states[5]
    Kpq = states[11]

    M = U/a # Mach number
    beta = sqrt(1 - M^2) #Prandtl-Glauert compressibility correction factor
    TI = c/a




    ##### Prepare inputs for normal force coefficient. 
    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2)  
    k_alpha = 1/bot #Equation 1.11a #Checked that uses slope. 
    Talpha = 3*k_alpha*TI/4 #Equation 1.10a
    Nnoncirc_a = 4*Talpha*(Ka-Kpa)/M #Equation 1.18, Noncirculatory component of normal force due to changes in alpha

    # if c== 1.419
    #     @show Talpha, Ka, Kpa, M #Ka and Kpa are quite large. 
    # end

    bot = (1-M) + dcndalpha*M*M*beta*(A1*b1 + A2*b2) 
    k_q = 1/bot #Equation 1.11b #Checked that uses slope. 
    Tq = 3*k_q*TI/4 #Equation 1.10b
    Nnoncirc_q = Tq*(Kq - Kpq)/M #Equation 1.19

    Nnoncirc_aq = Nnoncirc_a + Nnoncirc_q #Equation 1.17, Noncirculatory component of normal force via superposition.

    # if c== 1.419
    #     @show Cnc_nalpha, Cnc_nq #Cnc_nalpha is quite large. 
    # end

    alpha = states[1]
    q = states[3]

    X1 = states[6]
    X2 = states[7]
    X3 = states[8] #EQ 1.16a 
    X4 = states[9] #EQ 1.16a
    alphae = (alpha - alpha0) - X1 - X2 #EQ 1.14, Effective angle of attack

    dcndalpha_circ = dcndalpha/beta #EQ 1.12, Circulatory component of the normal force coefficient response to step change in alpha. #Checked that uses slope. 
    # Ncirc_q = dcndalpha_circ*(q - X3 - X4)/2 #EQ 1.16 
    Ncirc_q = dcndalpha_circ*q/2 - X3 -X4 #EQ 1.16 as Openfast v3.3.0 shows it. #Todo: This can't be right. It's asking for the circulatory normal force due to pitching, and it's using states X3 and X4.. I thought those states were for the chordwise force.... 

    Ncirc_aq = dcndalpha_circ*alphae #+ Cc_nq #EQ 1.13, Circulatory component of normal force via lumped approach.

    



    fpp = states[19]
    # println("coefficients fpp: ", fpp)
    fterm = (1 + 2*sqrt(fpp))/3

    # Cfsn = Cnc_naq + Cc_naq*(fterm)^2 #+ Cc_nq #EQ 1.38, Normal force coefficient after accounting for separated flow from TE #Todo: I'm not sure that this matches Gonzalez's modifcations. #Note: There is a typo in their equation. The first part of 1.53b says that it is Cfsn, and the second part is missing the Cc_nq term from 1.39, thus either they added an extra term in 1.39, or they missed a term in 1.53b. -> It didn't make much of a difference in the NREL 5MW verification case. -> I need to check OpenFAST and see what they do. 
    # Cfsn = Nnoncirc_aq + Ncirc_q + dcndalpha_circ*alphae*(fterm^2) #EQ 1.39
    Cfsn = Nnoncirc_aq + Ncirc_q + Ncirc_aq*(fterm^2) 

    # if c== 1.419
    #     @show Cnc_naq, Cc_naq, fterm #Both the Cnc and the Cc terms are large, the Cnc term is especially large. 
    # end

    Cvn = states[23]

    ### Total normal force 
    # Cn = Cfsn + Cvn #EQ 1.53b

    tauv = states[22]

    if tauv>0 #OpenFAST v3.3.0 - UnsteadyAero.f90 line 3230 #Note: Didn't appear to make any difference. 
        Cn = Cfsn + Cvn #EQ 1.53b
    else
        Cn = Cfsn
    end

    # if c== 1.419
    #     @show Cfsn, Cvn #Problem is Cfsn
    # end



    ######### Prepare inputs for chordwise force coefficient. 
    # Cpot = Ncirc_aq*tan(alphae + alpha0) #Equation 1.21   #Todo. What is Npot,circ? -> It's Ncirc_aq
    Cpot = dcndalpha_circ*alphae*alpha #What OpenFAST v3.3.0 has. #Todo. What is alpha? -> It appears to be the unfiltered angle of attack #Todo: What is C_nalpha_circ? -> At first I thought I had misread... but no... they definitely have C_nalpha_circ.... which suggests that maybe they mistyped. 

    fppc = states[21]
    # @show eta, Cpot
    Cfsc = Cpot*eta*(sqrt(fppc)-0.2) #EQ 1.40, Gonzalez modifications #Todo. This modification is pushing it crazy off. -> This was pushing it crazy off, but when I switched to using fppc it worked. 
    # Cfsc = Cpotc*eta*(sqrt(fppc))
    # Cfsc = Cpotc*eta*(sqrt(fpp))

    ### Chordwise force 
    Cc = Cfsc #EQ 1.55b




    ### Lift and Drag
    Cl = Cn*cos(alpha) + Cc*sin(alpha)
    Cd = Cn*sin(alpha) - Cc*cos(alpha) + Cd0 #Adding frictional drag back in. 



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
    Cn = Cl*cos(alpha) + Cd*sin(alpha) #Cd has friction added back in. 
    Cc = Cl*sin(alpha) - Cd*cos(alpha)
    Cm = 1

    return Cn, Cc, Cl, Cd, Cm
end

function initialize_ADG(Uvec, aoavec, tvec, airfoil::Airfoil, c, a) 
    # Cnfit = airfoil.cl
    dcndalpha = airfoil.dcndalpha
    alpha0 = airfoil.alpha0
    A1, A2 = airfoil.A
    b1, b2 = airfoil.b
    Tp, Tf0, Tv0, Tvl = airfoil.T
    _, alpha1 = airfoil.alphasep
    xcp = airfoil.xcp

    Cd0 = airfoil.cd(alpha0)

    Cn1 = airfoil.cl(alpha1)*cos(alpha1) + (airfoil.cd(alpha1) - Cd0)*sin(alpha1) #Removing the effect of frictional drag
    Cm0 = airfoil.cm(alpha0)



    dt = tvec[2] - tvec[1]

    aoa = aoavec[1]
    alpha = alphaf = aoavec[1] #No delay to begin. 
    aoadot = q = Ka = 0 #(aoavec[2] - aoavec[1])/dt #Ka was initialized too high. aoadot isn't used. I might be able to use q = (aoavec[2] - aoavec[1])*c/(U*deltat) and ka = q*U/deltat... Why are q and Ka both states?... if one is just a multiple of the other. I guess U and delta t change, so it isn't the same multiple across time. 

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

    states = [alpha, alphaf, q, Ka, Kq, X1, X2, X3, X4, Kpa, Kpq, Kppq, Kpppq, Dp, Df, Dfc, Cpotn, fp, fpp, fpc, fppc, tauv, Cvn, Cv, LESF, TESF, VRTX]

    Cn = airfoil.cn(alpha) 
    Cc = airfoil.cc(alpha)

    loads = [Cn, Cc, airfoil.cl(aoa), airfoil.cd(aoa), airfoil.cm(aoa)]

    # println("Df initialize: ", Df)


    p = [c, a, dcndalpha, alpha0, Cd0, Cm0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, xcp] #16 elements (-> total of 18 elements in p. )

    envvars = [Uvec[1], aoavec[1]]

    return states, loads, vcat(p, envvars)
end

function updateenvironment_ADG(p, U, aoa) #Todo:
    p[17] = U
    p[18] = aoa
end