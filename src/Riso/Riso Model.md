# Risø Model

The Risø model is a simplified version of the Beddoes-Leishman model, designed specifically for wind turbines. It has four states. The four states are:

|     Variable     |                             Name                             |                           Comments                           |
| :--------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
|      $x_1$       |                                                              | Using Theodorsen theory, we can come up with two ODEs to describe the change in the wake. |
|      $x_2$       |                                                              | The first state and second state model these changes. One is for high frequency changes while the other is for low frequency changes. |
| $x_3 = C_L^{p'}$ | Delayed coefficient of lift due to pressure. Delayed unsteady attached lift coefficient. | A simple time lag (first order lag) between the pressure field and the lift. |
|   $x_4 = f''$    |     Delayed seperation point? Unsteady seperation point.     |           I believe that this should range (0,1).            |







## State Rate Equations

$$
\begin{aligned}
\dot{x}_1 = - \frac{x_1}{T_u}\left(b_1 + \frac{c\dot{U}}{2U^2}\right) + \frac{b_1A_1 \alpha_{3/4}}{T_u} \\
\dot{x}_2 = - \frac{x_2}{T_u}\left(b_2 + \frac{c\dot{U}}{2U^2}\right) + \frac{b_2A_2 \alpha_{3/4}}{T_u} \\
\dot{x}_3 = - \frac{x_3}{T_p} + \frac{1}{T_p}\left(\frac{\partial C_l}{\partial \alpha}(\alpha_E - \alpha_0) + \pi T_u \dot{\alpha} \right) \\
\dot{x}_4 = - \frac{x_4}{T_f} = \frac{1}{T_f}f^{st}\left(\alpha_f\right)
\end{aligned}
$$



where $x_i$ are the states of the model, $T_x$ are time constants, $A_i$ and $b_i$ are airfoil specific constants describing the response of the airfoil, $U$ is the freestream velocity, $\alpha_{3/4}$ is the angle of attack the three-quarter point, $\frac{\partial C_l}{\partial \alpha}$ is the slope of the static lift curve in the linear region, $\alpha_E$ is the equivalent angle of attack, $f^{st}$ is the seperation point function, and $\alpha_f$ is the delayed angle of attack (this is something that I made up to make the notation a little clearer. ). 

A couple of notes: first, $T_u$ is calculated, whereas the other two time constants are specific for the airfoil (although several authors have said that they are pretty constant from airfoil to airfoil... but at the same time, those authors give different values for the time constant than the other authors... so it seems pretty arbitrary). Second, variables with a dot above it ($\dot{U}$ for instance) indicate the derivative with respect to time. Third, recognize that $U$ and $\alpha$ are both functions of time. Fourth, I just use the angle of attack for $\alpha_{3/4}$... which appears to work, but Hansen's paper specifically seperates from the geometric angle of attack.

Now I'm going to dive into a bunch of the functions. 

##### First Time Constant

$$
T_u = \frac{c}{2U}
$$

where $c$ is the chord, and $U$ is the freestream velocity. Since $U$ is a function of time, then this time constant isn't really constant...

##### Lift Curve Slope

Hansen gives an equation for $\frac{\partial C_l}{\partial \alpha}$, which makes it so that $f^{st}$ never goes over 1.... which is just a max equation... but I feel like it is a lame excuse and it is more accurate if you use a value that matches the majority of the linear region. I feel like there should be a better way to define $f^{st}$ so it never goes over 1, and plays nicely with optimization. At the same time he says: "However, in cases where lift curves are provided directly from measurements, or CFD computations, it can be necessary to use linear regression to determine a linear lift curve in a range of lower angles of attack with fully attached flow." DG has said several times that it is fairly common to choose constants that recreate the static lift curve. 
$$
\frac{\partial C_l}{\partial \alpha} = \text{max}\left(\frac{C_l^{st}(\alpha)}{(\alpha - \alpha_0)}\right)
$$
where $C_l^{st}$ is the lift from the static lift curve, and $\alpha_0$ is the zero lift angle of attack. You'll frequently see $\frac{\partial C_l}{\partial \alpha}\left(\alpha - \alpha_0\right)$, which is just the inviscid lift. Note that in his paper, Hansen denotes the linear lift curve slope as $C_{L,\alpha}$. 

##### Equivalent Angle of Attack

$$
\alpha_E = \alpha_{3/4}(1 - A_1 - A_2) + x_1 + x_2
$$

It's an equivalent angle of attack. Basically from what I understand is it is a weighting of what we expect the angle of attack to be from the geometric angle of attack and what the wake currently looks like.

##### Seperation Point Function

Here is a function that gives me a little bit of problem:
$$
f^{st}(\alpha) = \left(2 \sqrt{\frac{C_l^{st}}{C_l^i}} -1\right)^2
$$


where $C_l^i$ is the inviscid lift coefficient. The inviscid lift coefficent, as previously mentioned, is simply: $C_l^i(\alpha) = \frac{\partial C_l}{\partial \alpha} (\alpha - \alpha_0)$. However, the value of this should always be on the range (0,1). It should be equal to 1 when the flow is fully attached (in the linear region), and 0 when the angle of attack is outside the angles of full seperation. Hansen defines the angles of full seperation as when the static lift is one-fourth the value of the inviscid lift:
$$
|C_l^{st}(\alpha^{\pm fs})| = \left|\frac{\partial C_l}{\partial \alpha}\frac{(\alpha^{\pm fs} - \alpha_0)}{4}\right|
$$
What I currently have implemented for the seperation pooint function is: 
$$
f^{st}(\alpha)=\begin{cases}
          0 \quad &\text{if } \, \alpha<\alpha^{-fs} ~||~ \alpha>\alpha^{+fs} \\
          1 \quad &\text{if } \, f^{st} > 1 ~ || ~ f^{st} = \text{NaN} \\
          \left(2 \sqrt{\frac{C_l^{st}}{C_l^i}} -1\right)^2 \quad & \text{else}
     \end{cases}
$$
Which isn't my favorite. But it works. 



##### Other Time Constants

$T_f$ describes the time lag in the boundary layer. 



## Equations for the Coefficients

