
function Get_ExpEqn( dt, T, Y_minus1, x, x_minus1 )
# ! Called by : ComputeKelvinChain
# ! Calls  to : NONE
# !..............................................................................

   # dt       !< numerator of the exponent, \f$\Delta t \f$
   # T        !< denominator of the exponent (time constant), \f$T\f$
   # Y_minus1 !< previous value of the computed decay function, \f$ Y_{n-1} \f$
   # x        !< current value of x, \f$x_n\f$
   # x_minus1 !< previous value of x, \f$x_{n-1}\f$
   
   
   tmp = -dt/T # tmp should always be negative... should we check this, or are there some physics that make this check unnecessary?
   
   Get_ExpEqn = Y_minus1*exp(tmp) + (x - x_minus1)*exp(0.5*tmp)

end 


function ComputeKelvinChain( i, j, u, p, xd, OtherState, misc, AFInfo, KC, BL_p, ErrStat, ErrMsg )

    # dynamicFilterCutoffHz                         # find frequency based on reduced frequency of k = filtCutOff
       
    fprimeprime_m = 0
    
    M = U / a_s

    if M > 1 #call UA_CheckMachNumber(M, misc%FirstWarn_M, ErrStat2, ErrMsg2 ) #TODO: I haven't looked at the limit they set
        @warn("Mach number superseded 1.")
    end
          
    beta_M_Sqrd = 1.0 - M**2
    beta_M = sqrt(beta_M_Sqrd) 
       
    ##### Lookup values using Airfoil Info module
    #    call AFI_ComputeUACoefs( AFInfo, u%Re, u%UserProp, BL_p, ErrMsg2, ErrStat2 )
    #       call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
    #       if (ErrStat >= AbortErrLev) return
       
    C_nalpha_circ  =  C_nalpha / beta_M
          
    ## Override eta_e if we are using Flookup
    if FLookup
        eta_e = 1.0
    end
       
    ##### Compute Kelvin chain
    ds = 2.0*U*dt/c[i,j]                            # Eqn 1.5b
       
    if FirstPass[i,j]
       alpha_minus1 = alpha      
       alpha_filt_minus1 = u%alpha     
    else
       alpha_minus1 = alpha_minus1[i,j]
       alpha_filt_minus1 = alpha_filt_minus1[i,j]
    end
       
    #    # Low Pass filtering of alpha, q, and Kq in order to deal with numerical issues that arise for very small values of dt
    #    # See equations 1.7 and 1.8 of the manual
    #    # This filter is a Simple Infinite Impulse Response Filter
    #    # See https://en.wikipedia.org/wiki/Low-pass_filter#Simple_infinite_impulse_response_filter
       
    dynamicFilterCutoffHz = max( 1.0, U ) * filtCutOff / PI / c[i,j]
    LowPassConst  =  exp(-2.0*PI*dt*dynamicFilterCutoffHz) ## from Eqn 1.8 [7]
       
    alpha_filt_cur = LowPassConst*alpha_filt_minus1 + (1.0-LowPassConst)*alpha ## from eq 1.8 [1: typo in documentation, though]
       
    dalpha0  = alpha_filt_cur - alpha0
       
        
    ##### Compute Kalpha using Eqn 1.7
      
    Kalpha = (alpha_filt_cur - alpha_filt_minus1)/dt ## Eqn 1.7, using filtered values of alpha -> This is another route rather than using the equations presented in eq 1.8
       
    if FirstPass[i,j]
       Kalpha_f_minus1 = 0.0
    else
       Kalpha_f_minus1 = Kalpha_f_minus1[i,j]
    end
       
    q_cur = Kalpha*c[i,j] / U   #Kalpha here uses the low-pass filtered alphas (Eqn 1.8 [2])
       
           
    if FirstPass[i,j]
       q_minus1   = q_cur 
       q_f_minus1 = q_cur    
    else    
       q_minus1   = q_minus1[i,j]
       q_f_minus1 = q_f_minus1[i,j]
    end
       
    q_f_cur = LowPassConst*q_f_minus1 + (1.0-LowPassConst)*q_cur # (Eqn 1.8 [3])
 
    Kalpha_f = q_f_cur * U / c[i,j]  # Kalpha_f is using the low-pass filtered q (Eqn. 1.8 [4])
    
       
    #Compute Kq  using Eqn 1.8  with time-shifted q s
    #ifdef TEST_THEORY
    # Kq = (q_f_cur  - q_f_minus1 ) / dt #bjj: jmj questions if this should be the way it's implemented
    #else
       
    Kq = (q_cur  - q_minus1) / dt # Eqn. 1.8 [5]
    
    
    if FirstPass[i,j]
       Kq_f_minus1 = 0.0
    else
       Kq_f_minus1 = Kq_f_minus1[i,j]
    end
    
    Kq_f =  LowPassConst*Kq_f_minus1 + (1.0-LowPassConst)*Kq #Eqn. 1.8 [6]
       
       
    #bjj: todo: should we check the denominator to make sure it doesn't go to 0?
    
    k_alpha = 1.0 / ( (1.0 - M) + (C_nalpha/2.0) * M**2 * beta_M * (A1*b1 + A2*b2) )  # Eqn 1.11a
    k_q = 1.0 / ( (1.0 - M) +  C_nalpha * M**2 * beta_M * (A1*b1 + A2*b2) )  # Eqn 1.11b   
       T_I        = c(i,j) / a_s                                                                                                # Eqn 1.11c
       
    T_alpha  = T_I * k_alpha * 0.75    # Eqn 1.10a
    T_q      = T_I * k_q * 0.75 # Eqn 1.10b
          
    T_f           = T_f0 / sigma1[i,j]    # Eqn 1.37    
    T_fc          = T_f0 / sigma1c[i,j]    # NOTE: Added equations for time constants of fc (for Cc) and fm (for Cm) with UAMod=2
    T_fm          = T_f0 / sigma1m[i,j]
      
    Kprime_alpha  = Get_ExpEqn(dt, T_alpha, Kprime_alpha_minus1[i,j], Kalpha_f, Kalpha_f_minus1 )    # Eqn 1.18b

    Cn_alpha_nc   = 4.0*T_alpha * ( Kalpha_f - Kprime_alpha ) / M       # Eqn 1.18a
       
    Kprime_q = Get_ExpEqn(dt, T_q, Kprime_q_minus1[i,j], Kq_f, Kq_f_minus1)    # Eqn 1.19b 
    Cn_q_nc = -1.0*T_q * ( Kq_f - Kprime_q ) / M     # Eqn 1.19a
             
    Cn_alpha_q_nc = Cn_alpha_nc + Cn_q_nc                                                                             # Eqn 1.17
       
    if ShedEffect
       X1 = Get_ExpEqn( ds*beta_M_Sqrd*b1, 1.0, X1_minus1[i,j], A1*(alpha_filt_cur - alpha_filt_minus1), 0.0 ) # Eqn 1.15a
       X2 = Get_ExpEqn( ds*beta_M_Sqrd*b2, 1.0, X2_minus1[i,j], A2*(alpha_filt_cur - alpha_filt_minus1), 0.0 ) # Eqn 1.15b
    else
       X1  = 0.0  # u%alpha (and alpha_filt_cur) contains shed vorticity effect already 
       X2  = 0.0  # so that alpha_e = u%alpha-alpha0 directly
    end
       
    alpha_e = (alpha_filt_cur - alpha0) - X1 - X2   # Eqn 1.14
       
    Cn_alpha_q_circ = C_nalpha_circ * alpha_e   # Eqn 1.13
       
    if UAMod == UA_Gonzalez
             # Compute X3 and X4 using Eqn 1.16a  and then add Cn_q_circ (Eqn 1.16) to the previously computed Cn_alpha_q_circ
        if (ShedEffect) then
            X3              = Get_ExpEqn( ds*beta_M_Sqrd*b1, 1.0, X3_minus1(i,j), A1*(q_f_cur - q_f_minus1), 0.0 ) # Eqn 1.16a [1]
             X4              = Get_ExpEqn( ds*beta_M_Sqrd*b2, 1.0, X4_minus1(i,j), A2*(q_f_cur - q_f_minus1), 0.0 ) # Eqn 1.16a [2]
        else
            X3 = 0.0 # Similar to X1 and X2, we assumed that this effect is already included
            X4 = 0.0
        end
          
        Cn_q_circ       = C_nalpha_circ*q_f_cur/2.0 - X3 - X4                                                    # Eqn 1.16
    
    else # these aren't used (they are possibly output to UA output file (when UA_OUTS defined) file, though)
        X3              = 0.0
        X4              = 0.0
        Cn_q_circ       = 0.0
    end
       
    K3prime_q = Get_ExpEqn( b5*beta_M_Sqrd*ds, 1.0, K3prime_q_minus1[i,j],  A5*(q_f_cur - q_f_minus1), 0.0 )  # Eqn 1.26
    Cm_q_circ = -C_nalpha*(q_f_cur - K3prime_q)*c(i,j)/(16.0*beta_M*u%U)    # Eqn 1.25
       
    Cn_pot   = Cn_alpha_q_circ + Cn_alpha_q_nc   # Eqn 1.20 [2a]
       
       k_mq            = 7.0 / (15.0*(1.0-M) + 1.5 * C_nalpha * A5 * b5 * beta_M * M**2)       # Eqn 1.29 [2]      # CHECK THAT DENOM ISN'T ZERO#
       Kprimeprime_q   = Get_ExpEqn( real(dt,), k_mq**2*T_I , Kprimeprime_q_minus1(i,j) ,  Kq_f , Kq_f_minus1  )      # Eqn 1.29 [3]
       
          # Compute Cm_q_nc 
       if ( UAMod == UA_MinnemaPierce ) then
          Cm_q_nc =  -1.0 * Cn_q_nc / 4.0 - (k_alpha**2) * T_I * (Kq_f - Kprimeprime_q) / (3.0*M)      # Eqn 1.31
       else  
          Cm_q_nc = -7.0 * (k_mq**2) * T_I * (Kq_f - Kprimeprime_q) / (12.0*M)                                    # Eqn 1.29 [1]       
       end if
       
       if ( UAMod == UA_Gonzalez ) then
          Cc_pot = C_nalpha_circ * alpha_e * u%alpha  #Added this equation with (u%alpha) instead of tan(alpha_e+alpha0). First, tangent gives problems in idling conditions at angles of attack of 90 degrees. Second, the angle there is a physical concept according to the original BL model, and u%alpha could be more suitable 
       else   
       ### THIS IS A PROBLEM IF alpha_e+alpha0 ARE NEAR +/-PI/2
          Cc_pot = Cn_alpha_q_circ * tan(alpha_e+alpha0)                                                           # Eqn 1.21 with cn_pot_circ=Cn_alpha_q_circ as from Eqn 1.20 [3]  
       endif
       
       if (FirstPass(i,j)) then
          Cn_pot_minus1 = Cn_pot
       else
          Cn_pot_minus1 = Cn_pot_minus1(i,j)
       end if
       
       Dp            = Get_ExpEqn( ds, T_p, Dp_minus1(i,j), Cn_pot, Cn_pot_minus1 )                              # Eqn 1.35b
       Cn_prime      = Cn_Pot - Dp                                                                                       # Eqn 1.35a
       
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    # This code is taken from ADv14 but doesn't reflect the original intent of the UA theory document
    #ifdef TEST_THEORY
       if (FirstPass(i,j)) then
          Cn_prime_diff = 0.0
       else
          Cn_prime_diff = Cn_prime - Cn_prime_minus1(i,j)
       end if
          
    IF ( UAMod /= UA_Gonzalez ) THEN
       IF ( alpha_filt_cur * Cn_prime_diff < 0. ) THEN
    
          T_f   = T_f0*1.5
       ELSE
    
          T_f   = T_f0
       ENDIF
    ENDIF
    #endif   
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    
       alpha_f       = Cn_prime / C_nalpha_circ + alpha0                                                            # Eqn 1.34
       
       if (flookup) then
          call Get_f_from_Lookup( UAMod, u%Re, u%UserProp, alpha_f, alpha0, C_nalpha_circ, AFInfo, ErrStat2, ErrMsg2, f=fprime)     # Solve Eqn 1.32a for f (=fprime) when alpha is replaced with alpha_f (see issue when C_nalpha_circ is 0) 
          call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
          if (ErrStat >= AbortErrLev) return
       else   
          fprime = Get_f( alpha_f, alpha0, alpha1, alpha2, S1, S2, S3, S4)               # Eqn 1.33
       end if
       
       if (FirstPass(i,j)) then
          Df = 0.0
       else
          Df = Get_ExpEqn( ds, T_f, Df_minus1(i,j), fprime, fprime_minus1(i,j) )                                # Eqn 1.36b
       end if
          
       fprimeprime   = fprime - Df                                                                                       # Eqn 1.36a
       
       if (Flookup) then
             # Compute fprime using Eqn 1.32 and Eqn 1.33
          fprime_c   = Get_f_c_from_Lookup( UAMod, u%Re, u%UserProp, alpha_f, alpha0, C_nalpha_circ, eta_e, AFInfo, ErrStat2, ErrMsg2) # Solve Eqn 1.32b for f when alpha is replaced with alpha_f
             call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
             if (ErrStat >= AbortErrLev) return
    
          if ( UAMod == UA_Gonzalez ) then   #Added this part of the code to obtain fm
             call AFI_ComputeAirfoilCoefs( alpha_f, u%Re, u%UserProp, AFInfo, AFI_interp, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               if (ErrStat >= AbortErrLev) return
    
             Cn_temp = AFI_interCl*cos(alpha_f) + (AFI_interCd - AFI_interCd0)*sin(alpha_f)
             if (abs(Cn_temp) < 0.01 ) then
                fprime_m = 0.0
             else
                fprime_m = (AFI_interCm - AFI_interCm0) / Cn_temp 
             end if
          else
             fprime_m = 0.0
          endif
    
         
          if (FirstPass(i,j)) then
             Df_c = 0.0
             Df_m = 0.0
          else
             Df_c = Get_ExpEqn( ds, T_fc, Df_c_minus1(i,j), fprime_c, fprime_c_minus1(i,j)  )
             Df_m = Get_ExpEqn( ds, T_fm, Df_m_minus1(i,j), fprime_m, fprime_m_minus1(i,j)  )  # used in UAMod=UA_Gonzalez only
          end if
       
             # Compute Df using Eqn 1.36b   
       
             # Compute fprimeprime using Eqn 1.36a
          fprimeprime_c   = fprime_c - Df_c
          
          IF ( UAMod == UA_Gonzalez ) THEN
             fprimeprime_m   = fprime_m - Df_m
          END IF
       else
          fprime_c = fprime
          Df_c = Df
          fprimeprime_c = fprimeprime
    
             # variables used for UAMod=UA_Gonzalez
          fprime_m = 0.0
          Df_m     = 0.0
          fprimeprime_m = fprimeprime
       end if
       
       
       if ( UAMod == UA_Gonzalez ) then
          Cn_FS   = Cn_alpha_q_nc + Cn_q_circ + Cn_alpha_q_circ *  ( (1.0 + 2.0*sqrt(fprimeprime) ) / 3.0 )**2     # Eqn 1.39
       else
       # change proposed by Pariya:
         # Cn_FS   = Cn_alpha_q_nc                + Cn_alpha_q_circ *  ( (1.0 +          sqrt(fprimeprime) ) / 2.0 )**2     # Eqn 1.38
          call Get_f_from_Lookup( UAMod, u%Re, u%UserProp, alpha_e+alpha0, alpha0, C_nalpha_circ, AFInfo, ErrStat2, ErrMsg2, cn_fs=Cn_fs_temp)
          Cn_FS  = Cn_alpha_q_nc + C_nalpha_circ * alpha_e*fprimeprime  + Cn_fs_temp*(1-fprimeprime)
       end if
          
       if ( UAMod == UA_MinnemaPierce ) then
          if (FirstPass(i,j)) then     
             Dalphaf    = 0.0
          else
             Dalphaf    = Get_ExpEqn( ds, 0.1*T_f, Dalphaf_minus1(i,j), alpha_f, alphaf_minus1(i,j) )         # Eqn 1.43
          end if
       else
          Dalphaf    = 0.0
       end if
       
       if ( UAMod == UA_Gonzalez ) then
          C_V   = Cn_alpha_q_circ * ( 1.0 - ((1.0 + 2.0*sqrt(fprimeprime) )/3.0)**2  )               # Eqn. 1.50
       else
          C_V   = Cn_alpha_q_circ *  ( 1.0 - ( 0.5 + 0.5*sqrt(fprimeprime) )**2 )                         # Eqn. 1.49
       end if
    
       T_V = T_V0 / sigma3(i,j)                                                                                # Eqn 1.48
       
       if (FirstPass(i,j)) then
          Cn_v = 0.0
       else
          if (tau_V(i,j) > T_VL .AND. Kalpha_f * dalpha0 > 0 ) then # .AND. (.not. LESF)
             # We no longer require that T_V will always equal T_V0/2 when this condition is satisfied as was the case in AD v13 GJH 7/20/2017
             # If we fall into this condition, we need to require we stay here until the current vortex is shed (i.e., tauV is reset to zero)
             if ( UAMod == UA_Gonzalez ) then   #Added this equation from the formulation used in UAMod=UA_Gonzalez
                Cn_v = Cn_v_minus1(i,j)*exp(-2.0*ds/T_V)
             else    
                Cn_v = Cn_v_minus1(i,j)*exp(-ds/T_V)                                                                  # Eqn 1.52
             end if
          else      
             Cn_v = Get_ExpEqn( ds, T_V, Cn_v_minus1(i,j), C_V, C_V_minus1(i,j) )                               # Eqn 1.47
          end if
       
          if ( Cn_v < 0.0 ) then
             Cn_v = 0.0
          end if      
       end if
    
    end subroutine ComputeKelvinChain