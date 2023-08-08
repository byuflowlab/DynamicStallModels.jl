# Risø Theory

The Risø model$^2$ is a simplified version of the Beddoes-Leishman model, designed specifically for wind turbines.

It has four states. The four states are:

|     Variable     |                             Name                             |                           Comments                           |
| :--------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
|      $x_1$       | Effective Downwash Coefficient #1                                                             | Using Theodorsen theory, we can come up with two ODEs to describe the change in the wake. |
|      $x_2$       |   Effective Downwash Coefficient #2                                                            | The first state and second state model these changes. One is for high frequency changes while the other is for low frequency changes. |
| $x_3 = C_L^{p'}$ | Delayed coefficient of lift due to pressure. Delayed unsteady attached lift coefficient. | A simple time lag (first order lag) between the pressure field and the lift. |
|   $x_4 = f''$    |     Delayed seperation point movement.     |          This is caused by the dynamics of the boundary layer and should range (0,1).            |







## State Rate Equations

$$
\begin{aligned}
\dot{x}_1 = - \frac{x_1}{T_u}\left(b_1 + \frac{c\dot{U}}{2U^2}\right) + \frac{b_1A_1 \alpha_{3/4}}{T_u} \\
\dot{x}_2 = - \frac{x_2}{T_u}\left(b_2 + \frac{c\dot{U}}{2U^2}\right) + \frac{b_2A_2 \alpha_{3/4}}{T_u} \\
\dot{x}_3 = - \frac{x_3}{T_p} + \frac{1}{T_p}\left(\frac{\partial C_l}{\partial \alpha}(\alpha_E - \alpha_0) + \pi T_u \dot{\alpha} \right) \\
\dot{x}_4 = - \frac{x_4}{T_f} + \frac{1}{T_f}f^{st}\left(\alpha_f\right)
\end{aligned}
$$



where $x_i$ are the states of the model, $T_x$ are time constants, $A_i$ and $b_i$ are airfoil specific constants describing the response of the airfoil, $U$ is the freestream velocity, $\alpha_{3/4}$ is the angle of attack the three-quarter point, $\frac{\partial C_l}{\partial \alpha}$ is the slope of the static lift curve in the linear region, $\alpha_E$ is the effective angle of attack, $f^{st}$ is the seperation point function, and $\alpha_f$ is the delayed angle of attack.

A couple of notes: first, $T_u$ is calculated, whereas the other two time constants are specific for the airfoil. Second, variables with a dot above it ($\dot{U}$ for instance) indicate the derivative with respect to time. Third, recognize that $U$ and $\alpha$ are both functions of time. Fourth, DSM just uses the angle of attack for $\alpha_{3/4}$, which appears to work fine. However, Hansen's paper specifically seperates from the geometric angle of attack.

## Explanation of Variables

### First Time Constant

$$
T_u = \frac{c}{2U(t)}
$$

where $c$ is the chord length of the airfoil, and $U$ is the free stream velocity, which is a function of time.

### Second Time Constant

$$
T_p = \tau_p T_u
$$

where $T_u$ is the first time constant, and $\tau_p$ is a constant that is specific to the airfoil. In Faber's paper, he states that $\tau_p$ is defaulted to be at $1.5$, but this will generally change from case to case. This time constant is used to represent the pressure lag.

### Third Time Constant

$$
T_f = \tau_f T_u
$$

where $T_u$ is the first time constant, and $\tau_f$ has the same role as $\tau_p$ in the second time constant equation. However, Faber states that $\tau_f$ is defaulted to be $6.0$, but just as before, this value will generally vary from airfoil to airfoil. This constant is used to represent that lag in the boundary layer.

### Lift Curve Slope

Hansen gives an equation for $\frac{\partial C_l}{\partial \alpha}$, which makes it so that $f^{st}$ never goes over 1 (shown below). This is simply a max equation across the linear region for lift. Alternatively, a linear regression fit can be used across the inviscid region to determine a slope. If these options are not desirable, then simply using that $\frac{\partial C_l}{\partial \alpha} = 2\pi$ is also a possibility.

$$
\frac{\partial C_l}{\partial \alpha} = \text{max}\left(\frac{C_l^{st}(\alpha)}{(\alpha - \alpha_0)}\right)
$$

where $C_l^{st}$ is the lift from the static lift curve, and $\alpha_0$ is the zero lift angle of attack. You'll frequently see $\frac{\partial C_l}{\partial \alpha}\left(\alpha - \alpha_0\right)$, which is just the inviscid lift. Note that in his paper, Hansen denotes the linear lift curve slope as $C_{L,\alpha}$. 

### Equivalent Angle of Attack

$$
\alpha_E = \alpha_{3/4}(1 - A_1 - A_2) + x_1(t) + x_2(t)
$$

$\alpha_{3/4}$ is the angle of attack at the $3/4$ position on the chord, $A_1$ and $A_2$ are airfoil specific constants, and $x_1(t)$ and $x_2(t)$ are the first and second state values, respectively. The use for the effective angle of attack is that it combines the geometric contribution for angle of attack with pitching and plunging effects.

### Seperation Point Function


$$
f^{st}(\alpha) = \left(2 \sqrt{\frac{C_l^{st}}{C_l^i}} -1\right)^2
$$


where $C_l^i$ is the inviscid lift coefficient and $C_l^{st}$ is the static lift coefficient. The inviscid lift coefficent, as previously mentioned, is simply: $C_l^i(\alpha) = \frac{\partial C_l}{\partial \alpha} (\alpha - \alpha_0)$. The value of this should always be on the range (0,1). It should be equal to 1 when the flow is fully attached (in the linear region), and 0 when the angle of attack is outside the angles of full seperation. Hansen defines the angles of full seperation as when the static lift is one-fourth the value of the inviscid lift:

$$
|C_l^{st}(\alpha^{\pm fs})| = \left|\frac{\partial C_l}{\partial \alpha}\frac{(\alpha^{\pm fs} - \alpha_0)}{4}\right|
$$

what is implemented in the DSM code is the following piece wise function:

$$
f^{st}(\alpha)=\begin{cases}
          0 \quad &\text{if } \, \alpha<\alpha^{-fs} ~||~ \alpha>\alpha^{+fs} \\
          1 \quad &\text{if } \, f^{st} > 1 ~ || ~ f^{st} = \text{NaN} \\
          \left(2 \sqrt{\frac{C_l^{st}}{C_l^i}} -1\right)^2 \quad & \text{else}
     \end{cases}
$$

where a root finding function is used to find the angles of attack that satisfy the equation for the fully separated angle of attacks.

### Delayed Angle of Attack

$$
\alpha_f = \frac{x_3}{C_{L,\alpha}} +\alpha_0
$$

this is the delayed angle of attack value that is used the fourth state rate equation. It has been denoted like this for simplicity.

## Fully Separated Lift

$$
C_L^{fs} = \frac{C_L^{st} - C_{L,\alpha}(\alpha - \alpha_0)f^{st}}{1 - f^{st}}
$$

This equation gives the fully separated lift coefficient, which is then used in finding the dyanamic lift coefficient. This coefficient represents what the lift coefficient would be in an environment where the airflow on the airfoil is completely detatched. It appears that when using this equation a problem will arise when the separation point approaches a value of $1$. However, in Hansen's paper he assures that as long as the equation for the separation point is substituted into the fully separated lift equation, then:

$$
\lim_{f^{st} \to 1}C_L^{fs}(\alpha)  = \frac{C_L^{st}(\alpha)}{2} ~~ \mathrm{and} ~~ \lim_{f^{st} \to 1}C_{L,\alpha}(\alpha - \alpha_0)  = C_L^{st}(\alpha)
$$

## Equations for the Coefficients

### Dynamic Lift

$$
C_L^{dyn} = C_{L,\alpha} (\alpha_E - \alpha_0)x_4 + C_L^{fs}(\alpha_E)(1-x_4) + \pi T_u \dot\alpha
$$

where $\dot\alpha$ is the rate of change of the angle of attack with respect to time. Note that many of these varibales in this equation are functions with respect to time.

### Dynamic Drag

$$
C_D^{dyn} = C_D^{st}(\alpha_E) + C_L^{dyn}(\alpha - \alpha_E) + (C_D^{st}(\alpha_E) - C_{D0})(\frac{\sqrt{f^{st}(\alpha_E)} - \sqrt{x_4}}{2} - \frac{f^{st}(\alpha_E) - x_4}{4})
$$

where $C_D^{st}$ is the static drag coefficient and $C_{D0}$ is the coefficient of drag at the zero lift angle of attack.

### Dynamic Moment

$$
C_M^{dyn} = C_M^{st}(\alpha_E) + C_L^{dyn}(a^{st}(x_4) - a^{st}(f^{st}(\alpha_E))) - \frac{\pi}{2} T_u \dot\alpha
$$

where $C_M^{st}$ is the static moment coefficient and $a^{st}$ can be represented by the equation below:

$$
a^{st}(x) = \frac{C_M^{st}(x) - C_{M0}}{C_L^{st}(x)}
$$

where $C_{M0}$ is the coefficient of moment at the zero lift angle of attack. This equation gives the relationship of how far from the quarter chord position the equivalent pressure center has shifted.
