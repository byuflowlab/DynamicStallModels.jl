# Beddoes-Leishman Theory



## State Rate Equations 
#### Equations

==TODO: Enumeration not working properly.===
1) reduced linear lift coefficient - high frequency (or low, don't know which) (EQ 13)

```math
\dot{c}_1 + \omega_1 c_1(t) = A_1 \dot{C}_{L0}(t)
```



2) reduced linear lift coefficient - low frequency (EQ 13)

- reduced linear lift coefficient - low frequency (EQ 13)

```math
\dot{c}_2 + \omega_2 c_2(t) = A_2 \dot{C}_{L0}(t)
```



- impulsive contributions diminishing in time due to wave propagation. (EQ A.2)

```math
\dot{c}_3 + \omega_5 c_3(t) = \frac{4}{M}A_3\dot{\alpha}
```

- ''   ''   ''   ''   ''   ''   ''  (EQ A.2)

```math
\dot{c}_4 + \omega_6 c_4(t) = \frac{1}{M}A_4\frac{c}{V}\ddot{\alpha}
```

- The retarded linear lift $C'_{L0,d}(t)$​​  is introduced as a delayed state variable of the linear lift $C_{L0,d}(t)$​​​, which gives a one-to-one correspondence between the pressure coefficient and the dynamic lift at changing pitch rates (EQ A.3)

```math
\dot{C}'_{L0,d}(t) = -\omega_7\bigg(C'_{L0,d}(t) - C_{L0,d}(t)\bigg)
```

Which when you plug everything in, then you get a  different equation: (EQ A.5)
```math
\dot{C}'_{L0,d}(t) = \omega_7\bigg(-c_1(t) - c_2(t) + c_3(t) + c_4(t) - C'_{L0,d}(t) + C_{L0}(\alpha)\bigg)
```


- dynamic attachment degree (EQ 16)

```math
\dot{f}_d = - \omega_3 \bigg(f_d(t)-f(\alpha)\bigg)
```



- Leading edge vortex contribution (EQ 22)

```math
\dot{C}_{L,v}(t) + \omega_4 C_{L,v}(t) = 
\begin{cases}
\Delta \dot{C}_L(t) \quad \textrm{for} \quad \tau < 1 \hspace{3pt} \textrm{and} \hspace{3pt} \dot{\alpha}>0 \\
0 \hspace{38pt} \textrm{otherwise}
\end{cases}
```



- Time Constant** (EQ 21)

```math
\dot{\tau} = 
\begin{cases}
\frac{V}{3c} \quad \textrm{for} \quad \alpha>\alpha_v,\\
0 \quad \textrm{otherwise}
\end{cases}
```

 ** This state variable was not included in the original formulation, however, they included a time constant $\tau$ that was the position of the leading edge vortex, and it moves with a speed of $\dot{\tau}$​​, so I figured, I'd make it a state variable that was really simple. Let the solver manage that information. I don't think that I could write something to manage it better. 




## Indicial Formulation (Original)



## Indicial Formulation (AeroDyn)$^3$

AeroDyn has it's own implementation of the Beddoes-Leishman model, with several variations here and there. The method is presented here. Frequently the method is presented in a way that makes sense theoretically, here however, we will present the method in the order that the functions are calculated. Note that the report that this is implementation is patterned off of gives a slightly different order for the functions. They are ammended here, considering that not all of the originally given equations were required to update the states.



#### States
AeroDyn lists a set of discrete states that are different from what a demarcated as the states from the original Beddoes-Leishman model. These states make sense when you consider simulation, as you need this information from step to step. I'll note that I had to add a couple of states to the original documentation, specifically $\tau_v$, $\sigma_1$, and $\sigma_3$. 

|     Variable     | Code Variable |  Name     |          Comments               |
| :--------------: | :----: |:----------------------------------------------------------: | :----------------------------------------------------------: |
|      $\alpha$       | `alpha` |Angle of attack     | This state is low-pass-filtered. |
| $\alpha_{f}$ | `alphaf` |Delayed effective Angle of Incidence (AOI) |  |
|   $q$    | `q`  | Pitch rate   |    This state is low-pass-filtered.     |
| $K_{\alpha}$ | `Ka` | $K_{\alpha}$ | This state is low-pass-filtered. |
| $K_{q_{lp}}$ | `Kq` | $K_q$ | This state is low-pass-filtered. | 
|  $X_1$  | `X1`  | Circulatory portion of the angle of attack |         |
|   $X_2$  | `X2` | Second circulatory portion of the angle of attack |        |
|   $K'_\alpha$   | `Kpa` | Deficiency function for noncirulatory component of normal force based on angle of attack  |      |
|   $K'_q$   | `Kpq` |Deficiency function for noncirulatory component of normal force based on pitching rate  |      |
|   $K''_q$   | `Kppq` | Deficiency function for noncirulatory component of moment  |      |
|   $K'''_q$   | `Kpppq` | Deficiency function for cirulatory component of moment  |      |
|   $D_p$   | `Dp` | Deficiency function for cirulatory component of normal force  |      |
|   $D_f$   | `Df` | Deficiency function for the separation point   |      |
|   $C_n^{pot}$   | `Cpotn`|  Total normal force under attached conditions |      |
|   $f'$   | `fp` | Effective seperation point |      |
|   $f''$   | `fpp` | Delayed effective seperation point |      |
|   $\tau_v$   | `tauv` | Nondimensional time that the leading edge vortex has been travelling. |      |
|   $C_n^v$   |  `Cvn` | Normal force coefficient due to accumulated vorticity |      |
|   $C_v$   | `Cv` | Coefficient of vorticity?  |      |
|   $\sigma_1$   | `sigma1` | Constant to modify $T_f$. |  Changes discontinuously.  |
|   $\sigma_3$   | `sigma3` | Constant to modify $T_v$. |  Changes discontinuously.  |

|     Variable     |                             Name                             |                           Comments                           |
| :--------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
|      $\alpha$       |  Angle of Attack  (AOA)  | I don't see where this is used (in the previous time step). |
|      $\alpha_{lp}$       |  Low-pass-filtered AOA      |  |
| $\alpha_{f}$ | Delayed effective Angle of Incidence (AOI) | I don't see where this is used (in the previous time step). |
|   $q$    |    Pitch rate   |      I don't see where this is used (in the previous time step).     |
|  $q_{lp}$ | low-pass-filtered $q$ | |
| $K_{\alpha_{lp}}$ | low-pass-filtered $K_{\alpha}$ | 
| $K_{q_{lp}}$ | low-pass-filtered $K_q$ | 
|  $X_1$  |     |         |
|   $X_2$  |     |        |
|   $X_3$   |     |      |
|   $X_4$   |     |      |
|   $K'_\alpha$   |   Deficiency function for noncirulatory component of normal force based on angle of attack  |      |
|   $K'_q$   |  Deficiency function for noncirulatory component of normal force based on pitching rate  |      |
|   $K''_q$   |   Deficiency function for noncirulatory component of moment  |      |
|   $K'''_q$   |   Deficiency function for cirulatory component of moment  |      |
|   $D_p$   |   Deficiency function for cirulatory component of normal force  |      |
|   $D_f$   |  Deficiency function for the separation point   |      |
|   $D_{f_c}$   |     |      |
|   $C_n^{pot}$   |     |      |
|   $f'$   |     |      |
|   $f'_c$   |     |      |
|   $f''$   |     |      |
|   $f''_c$   |     |      |
|   $\tau_v$   |     |      |
|   $C_n^v$   |     |      |
|   $C_v$   |     |      |
|   $D_{\alpha f}$   |     |      |




#### Logical Flags
```julia
if UAmod == 1
    # Closest model to the original Leishman-Beddoes formulation
elseif UAmod == 2
    # Modifications to the original model and simplifications following Gonzalex (2014)
        #= Equations
        1.58 -> 1.60
        1.56 -> 1.55b
        1.53 -> 1.53b
        1.49 -> 1.50
        1.41 -> 1.45
        1.38 -> 1.39
        add 1.16 to Ccnalphaq(s,M) Eq(1.13)
        =#
else #UAmod == 3
    # Modifications to the original model and simplifications following Pierce (1996) and Minnema (1998)
        #= Equations
        1.56 -> 1.55
        1.58 -> 1.59
        1.41 -> 1.43-1.44
        Modify 1.30 with 1.31
        =#
end

if flookup # == true
    # EQ 1.33 gets replaced by lookup values for f'n and f'c. Note that if UAmod == 2 or 3, the flag is automatically set to true. 
end
```

#### Algorithm to Update Discrete States
The AeroDyn documentation provides this algorithm (with some minor modifications included). 

1) Equation 1.11c - Some time constant
    - $T_I = \frac{c}{a}$
2) Equation 1.5b - Nondimensional distance
    - $\Delta s = \frac{2}{c}U(t)\Delta t$

- 1.7 - 1.8
    - (1.7) Nondimensional pitch rate
        - $q = \frac{\dot{\alpha}c}{U} \approx \frac{K_{\alpha n}c}{U}$ 
        - $K_{\alpha n} = \frac{ \alpha_n - \alpha_{n-1}}{\Delta t}$

    - (1.8) Applying a low pass filter to $K_\alpha$.
        - $\alpha_{lp_n} = C_{lp}\alpha_{lp_{n-1}} + (1 - C_lp)\alpha_n$ (low-pass-filtered $\alpha$)
        - $q_n = \frac{(\alpha_{lp_n} - \alpha_{lp_{n-1}})c}{U_n \Delta t}$ 
        - $q_{lp_n} = C_{lp}q_{lp_{n-1}} + (1 - C_{lp})q_n$ (low-pass-filtered $q$)
        - $K_{\alpha lp_n} = \frac{q_{lp_n}U_n}{\Delta t}$ (modified value of $K_\alpha$)
        - $K_{q_n} = \frac{q_n - q_{n-1}}{\Delta t}$ (backward finite difference of $q$)
        - $K_{qlp_n} = C_{lp}K_{q lp_{n-1}} + (1 - C_{lp})K_{q_n}$ (low-pass-filtered $k_q$)
        - $C_{lp} = e^{-2\pi \Delta t \zeta_{lp}}$ (low-pass-filter constant)
        - $\zeta_{lp}$ typically - 3 dB (low-pass-frequency cutoff)

**From here on the LP subscript is dropped with the understanding that quantities such as $\alpha$, $K_\alpha$, $q$, and $K_q$ denote filtered quantities**. 

- 1.11a (calculated solely at first time step)
    - $k_\alpha(M) = \frac{1}{(1 - M) + C_{n\alpha}(M)M^2\beta_M(A_1b_1 +A_2b_2)}$
    - Note that this is different from $K_\alpha$. It is a lowercase k. 
    - Suggested to be calculated solely at first time step.
5) Equation 1.11b 
    - $k_q(M) = \frac{1}{(1-M) + C_{n\alpha}(M)M^2\beta_M(A_1b_1 + A_2b_2)}$
    - ==There appears to be a repeat of the equation.==
        - In Leishman's 1990 paper, the second term in the denominator of $k_\alpha$ is half of the $k_q$ term. 

- 1.10
    - $T_\alpha(M) = 0.75 k_\alpha(M)T_I$
    - $T_q(M) = 0.75 k_q(M)T_I$
7) Equation 1.37
    - $T_f = \frac{T_{f0}}{\sigma_1}$
8) Equation 1.48
    - $T_v = \frac{T_{v0}}{\sigma_3}$
9) Equation 1.18 - Noncirculatory component of normal force due to changes in $\alpha$
    - $C_{n_\alpha}^{nc}(s,M) = \frac{4T_\alpha(M)}{M}(K_\alpha - K'_\alpha)$
    - (Deficiency function for $C_{n_\alpha}^{nc}(s,M)$) 
        - $K'_{\alpha_n} = K'_{\alpha_{n-1}} \text{exp}\left(-\frac{\Delta t}{T_\alpha(M)} \right) + (K_{\alpha_n} - K_{\alpha_{n-1}})\text{exp}\left(- \frac{\Delta t}{2T_\alpha(M)} \right)$
10) Equation 1.19 - Noncirculatory component (of Normal force?) due to changes in $q$
    - $C_{n_q}^{nc}(s,M) = \frac{T_q(M)}{M}(K_{q_n} - K'_{q_n})$
        - Note that in the above equation a negative has been neglected. In the equation 1.20 the original documentation appeared to have dropped a negative that appears from this term and the solution comes out correct. 
    - $K'_{q_n} = K'_{q_{n-1}}\text{exp}\left(- \frac{\Delta t}{T_q(M)} \right) + (K_{q_n} - K_{q_{n-1}})\text{exp}\left(- \frac{\Delta t}{2 T_q(M)} \right)$
11) Equation 1.17 - Noncirculatory component of normal force via superposition)
    - $C_{n_{\alpha q}}^{nc}(s, M) = C_{n_\alpha}^{nc}(s,M) + C_{n_q}^{nc}(s,M)$
        - Note that the equation is stating that these normal force components are functions of $s$ and $M$. 
12) Equation 1.15 - Update States 1 and 2
    - $X_{1_n} = X_{1_{n-1}}\text{exp}(-b_1\beta_M^2\Delta s) + A_1 \text{exp}(-b_1\beta_M^2\Delta s/2)\Delta \alpha_n $
    - $X_{2_n} = X_{2_{n-1}}\text{exp}(-b_2\beta_M^2\Delta s) + A_2 \text{exp}(-b_2\beta_M^2\Delta s/2)\Delta \alpha_n $
13) Equation 1.14 - Effective angle of attack
    - $\alpha_{e_n}(s, M) = (\alpha_n - \alpha_0) - X_{1_n}(\Delta s) - X_{2_n}(\Delta S)$
        - I'm assuming that the equation is saying that the two states are functions of $\Delta s$. 

14) Equation 1.13 - Circulatory component normal force by the lumped approach
    - $C_{n_{\alpha, q}}^C(s, M) = C_{n\alpha}^C(s,M)\alpha_e$
        - Equation 1.12 $C_{n\alpha}^C(s,M) = \frac{C_{n\alpha}(M)}{\beta_M}$ 
15) Equation 1.26 
    - $K'''_{q_n} = K'''_{q_{n-1}}\text{exp}\left(-b_5 \beta_M^2 \Delta s\right) + A_5\Delta q_n \text{exp}\left(- b_5 \beta_M^2 \Delta s/2 \right)$
    - Note that $A_5 = 1$ and $b_5 = 5$.
16) Equation 1.20 - Total normal force under attached conditions
    - $C_n^{pot} = C_{n_{\alpha,q}}^C(s,M) + C_{n_{\alpha,q}}^{nc}(s,M)$
17) Equation 1.29 - Other noncirculatory component
    - $C_{m_q}^{nc}(s,M) = -\frac{7T_i}{12M}\left(k_{m,q}(M)\right)^2 (K_q - K''_q)$
        - $k_{m,q}(M) = \frac{7}{15(1-M) + 1.5 C_{n\alpha}(M)A_5 b_5 \beta_M M^2}$
        - $K''_{q_n} = K''_{q_{n-1}}\text{exp}\left(- \frac{\Delta t}{(k_{m,q}(M))^2 T_I} \right) + (K_{q_n} - K_{q_{n-1}})\text{exp}\left(- \frac{\Delta t}{2(k_{m,q}(M))^2 T_I} \right)$
            - ==Why are we calculating these portions of the moment?==

-  (Chordwise force)
    - (1.20c) $C^{pot,c}_n = C^c_{n\alpha}(s,M)\alpha_e$ 
        - Used above in equation 1.20 -> Could combine. 
    - (1.21) $C_C^{pot} = C_n^{pot,c}\text{tan}(\alpha_e + \alpha_0)$

- 1.35 (Lagged Circulatory normal force (boundary layer response))
    - $C'_n = C_n^{pot} - D_p$
        - $D_{p_n} = D_{p_{n-1}}\text{exp}\left(- \frac{\Delta s}{T_p}\right) + (C_{n_{n}}^{pot}-C_{n_{n-1}}^{pot})\text{exp}\left(-\frac{\Delta s}{2T_p}\right) $

- 1.34 (Delayed effective angle of incidence)
    - $\alpha_f = \frac{C'_n}{C_{n\alpha}^c(s,M)} + \alpha_0$

20) Equation 1.33 - TE Separation point distance (from LE in percentage chord)
```math
 f = \begin{cases}
    1 - 0.3 \text{exp}\left(\frac{\alpha - \alpha_1}{S_1}\right) \text{ if } \alpha_0 \leq \alpha \leq \alpha_1 \\
    1 - 0.3 \text{exp}\left(\frac{\alpha_2 - \alpha}{S_3}\right) \text{ if } \alpha_2 \leq \alpha \leq \alpha_0 \\
    0.04 + 0.66\text{exp}\left(\frac{\alpha_1 - \alpha}{S_2}\right) \text{ if } \alpha > \alpha _1 \\
    0.04 + 0.66 \text{exp}\left(\frac{\alpha - \alpha_2}{S_4}\right) \text{ if } \alpha < \alpha_2
    \end{cases}
```
- 
    - $S_1$ and $S_2$ are best-fit constants that define the abruptness of the static stall. 
    - $\alpha_1$ is the angle of attack @ $f=0.7$ for $\alpha\geq\alpha_0$. 
    - $\alpha_2$ is the angle of attack @ $f=0.7$ for $\alpha<\alpha_0$. 

21) Equation 1.36 - f that accounts for delays in the boundary layer
    - $f'' = f' - D_f$
        - $f' = f(\alpha_f)$
            - $f'$ can be derived from a direct lookup table of static airfoil data reversing EQ 1.32. In fact, two values of $f'$ could be calculated: one for $C_n(f'_n)$ and one for $C_c(f'_c)$. 
        - $D_{f_n} = D_{f_{n-1}} \text{exp}\left(- \frac{\Delta s}{T_f}\right) + (f'_n - f'_{n-1})\text{exp}\left(- \frac{\Delta s}{2T_f}\right)$

22) Equation 1.38 [or 1.39] - Normal force coefficient after accounting for separated flow from TE
    - (1.38) $C_n^{fs} = C_{n_{\alpha,q}}^{nc}(s, M) + C_{n_{\alpha,q}}^{c}(s,M)\left(\frac{1 + \sqrt{f''}}{2}\right)^2$
    - (1.39)
        - Gonzalez (2014) correction. 

23) Equation 1.49 [or 1.50] - Normal force coefficient due to accumulated vorticity
    - (1.49) $C_V = C_{n\alpha}^c(s,M)\alpha_e\left(1 - \frac{1 + \sqrt{f''}}{2}\right)^2$
    - (1.50)
        - Gonzalez (2014) correction.

24) Equation 1.47 [or 1.52] - Normal force coefficient contribution from the additional lift from LE vortex
    - (1.47) $C_{n_n}^v = C_{n_{n-1}}^v \text{exp}\left(- \frac{\Delta s}{T_v}\right) + (C_{V_n} - C_{V_{n-1}})\text{exp}\left(- \frac{\Delta s}{2T_v}\right)$
        - Note that $C_n^v$ is not allowed to have a sign opposite to that of $C_n^{fs}$. 
    - (1.52)
        - some correction

#### Update Other States
###### Leading Edge Separation
- if $C'_n > C_{n1} \implies (C'_n<C_{n2}$ for $\alpha<\alpha_0)$ 
    - `LESF = true #LE separation can occur`
- else
    - `LESF = false #Reattacment can occur`
- end
###### Trailing Edge Separation
- if $f''_t < f''_{t-1}$
    - `TESF = true #TE separation in progress `
- else
    - `TESF = false #TE reattachment in progress`
- end

###### Vortex advection
- if $0<\tau_v\leq 2T_{vl}$
    - `VRTX = true #Vortex advection in progress`
- else
    - `VRTX = false #Vortex is in wake`
- end

###### Vortex Position Reset
- if $\tau_v \geq 1 + \frac{T_{sh}}{T_{vl}}$ and `LESF = true`
    - $\tau_v=0$
- end

###### $T_f$ Modifications
eq 1.37
```math 
T_f = T_{f0}/\sigma_1 
```

$\sigma_1$ is initilized at 1. 
$\Delta_{\alpha 0} = \alpha - \alpha_0$

```julia
if TESF == true #(Separation)
    if K_alpha Delta_alpha0 < 0
        sigma_1 = 2 #(Accelerate separation point movement)
    else
        if LESF == false
            sigma_1 = 1 #(LE separation can occur)
        else
            if fpp_nm1 <= 0.7 
                sigma_1 = 2 #(accelerate separation point movement if separation is occuring )
            else
                sigma_1 = 1.75
            end
        end
    end
else #(reattachment (`TESF==false`))
    if LESF == false
        sigma_1 = 0.5 #(Slow down reattachment)
    elseif VRTX == true  && 0 <= tau_v <= Tvl
        sigma_1 = 0.25 #(No flow reattachment if vortex shedding is in progress)
    elseif K_alpha*Delta_alpha > 0
        sigma_1 = 0.75
    end 
end
```
Apparently Aerodyn uses a "simpler" version. 

###### $T_v$ Modifications

```julia
simga3 = 1

if Tvl <= Tau_v <= 2*Tvl
    sigma3 = 3 #Postshedding
    if TESF== false
        sigma3 = 4 #Accelerate vortex lift decay
        if VRTX == true && 0 <= Tau_v <= Tvl
            if K_alpha*Delta_alpha0 < 0
                sigma3 = 2 #Accelerate vortex lift decay
            else
                sigma3 = 1 #default
            end
        end
    end
else
    if K_alpha*Delta_alpha0 < 0
        sigma3 = 4 #vortex lift must decay fast
    end
end

if TESF == false && K_q*Delta_alpha0 < 0 #Note that it's Kq and not K_alpha
    sigma3 = 1 #Default
end

```

#### Calculate Output
==For the first time step, outputs are determined by static lookup tables==

**Same as Update Discrete States Plus**
- 1.53 (Total normal force)
    - $C_n = C_n^{fs} + C_n^v$
        - Does this somehow include a switch between separated flow and attached flow? Or do the components changed based on that? 
- 1.55 (Chordwise Force)
    - Equation 1.20c - $C^{pot,c}_n = C^c_{n\alpha}(s,M)\alpha_e$ 
        - Used above in equation 1.20 -> Could combine. 
    - Equation 1.21 - $C_C^{pot} = C_n^{pot,c}\text{tan}(\alpha_e + \alpha_0)$
    - $C_c = C_c^{fs} + C_n^v \text{tan}(\alpha_e)\left(1 - \frac{\tau_v}{T_{vl}}\right)$
        - In the current release $\text{tan}(\alpha_e)\approx \alpha_e$
- 1.2 (Lift and Drag)
    - $C_l = C_n \cos(\alpha) + C_c\sin(\alpha)$
    - $C_d = Cn\sin(\alpha) - C_c\cos(\alpha) + C_{d_0}$
- 1.58 or 1.59 (Minnema) or 1.60 (Gonzalez) (pitching moment)
    - Equation 1.25 - Circulatory component of moment
        - $C_{m_q}^c(s, M) = - \frac{C_{n\alpha}(M)}{16\beta_M} (q - K'''_q)\frac{c}{U}$
    - $C_{m_q}^{nc}(s,M) = -\frac{7T_i}{12M}\left(k_{m,q}(M)\right)^2 (K_q - K''_q)$
    - (1.58) $C_m = C_{m0} - C_{n_{\alpha,q}}^c(s,M)(\hat{x}_{cp} - 0.25) + C_{m_q}^c(s,M) + C_{m_\alpha}^{nc}(s,M) + C_{m_q}^{nc}(s,M) + C_m^v$
    - (1.59)
    - (1.60)