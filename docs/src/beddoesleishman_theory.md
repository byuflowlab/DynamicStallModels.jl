# Beddoes-Leishman Theory



## State Rate Equations 
###### Under Construction
#### Equations


-  reduced linear lift coefficient - high frequency (or low, don't know which) (EQ 13)

```math
\dot{c}_1 + \omega_1 c_1(t) = A_1 \dot{C}_{L0}(t)
```




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
###### Under Construction

## Indicial Formulation (AeroDyn)$^3$
###### Under construction
AeroDyn has it's own implementation of the Beddoes-Leishman model, with several variations here and there. I will not go through the theory here, as you can reference the original document. Note that the code is slightly different from the theory documentation. Some differences between the theory document and software include:
- Several calculations are used that are claimed to be more stable than the original theory (see the OpenFAST code).
- There are separate separation point and associated deficiency functions for the normal, chordwise, and moment loads. 
- There are values stored with the states that aren't labeled as states.

What I have implemented isn't quite the same as what AeroDyn presents in their theory document, but it produces the same loads (<0.01%, <0.05% relative error for the normal, tangential, and moment loads respectively) as OpenFAST v3.3.0. 

Some of the differences between my code and OpenFAST v3.3.0 include:
- I don't differentiate between other states and the continuous states. 
- I omit several states that don't need to be states. 
- I reorganize the states for readability. 
- There are some simple unstated values that I assume. 

If you want to see the exact algorithm that I implement, check out the code. 


AeroDyn also has several different versions that you can change between. The theory is also presented in the theory document mentioned above. Note that states and other variables change with the modifications (as expected). The different versions I have are:
- Original (Under construction)
- Gonzalez (Validated)






