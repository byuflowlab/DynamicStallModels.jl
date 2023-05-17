#  Øye Theory

Øye's model is a linear interpolation of fully attached lift (inviscid), and fully detached lift based on a dynamic degree of attachment. The dynamic attachment degree is a simple first-order filter on the separation point on the suction side of an airfoil. This dynamic degree of attachment is the Øye model's only state, which makes it the simplest dynamic stall model to implement.

The state is:

| Variable | Name | Comments |
| -------- | -------------------------- | ----- |
| $f^d$ | Dynamic Attachment Degree | This differential equation is also utilized in the Riso model and Beddoes-Leishman model |

## State Rate Equation

$$
\begin{equation}
\dot{f}^d = -\tau (f^d- f(\alpha))
\end{equation}
$$

As stated, equation 1 is the only state rate equation in this model. $f^d$ is the state, $\tau$ is the time constant, and $f^s(\alpha)$ is the position of the fluid's separation point as a function of the angle of attack (also known as the static separation point).

Additionally, it should be noted that variables with a dot above them, such as $\dot{f}^d$, are derivatives with respect to time.

The above equation can be evaluated using a prefered ODE solver. Once the solution is discovered, it can be implemented into the following equation to solve for the dynamic lift.

$$
\begin{equation}
C_{L}(t) = f^{d}(t)C_l^{inv} + (1-f^{d}(t))C_l^{fs}
\end{equation}
$$

The indicial approach simply solves equation 1 by assuming that the inflow velocity and angle of attack do not vary within the time step, and using an integrating factor:

$$
\begin{equation}
f^{d}_{i+1} = f^{s}(\alpha_i) + (f^{d}_i - f^{s}(\alpha_i))\hspace{3pt}e^{-\frac{dt}{\tau}}
\end{equation}
$$


## Explanation of Variables
#### Time Constant

The time constant is experimentally determined, and provides the needed information for how long it takes to come to steady state. 

$$
\begin{equation}
\tau = \frac{V}{A c}
\end{equation}
$$

where $\tau$ is the time constant, $V$ is the inflow velocity (m/s) $A$ is the airfoil time coefficient, and $c$ is the airfoil chord (m). $A$ is assumed to be $A=4$.
For myriad reasons, different authors have decided to rearrange this equation to suit their circumstance. Since we verified our code against Faber and Larsen, we present here their formulations for the airfoil time coefficient.

Faber multiplies the airfoil time coefficient by two:

$$
\begin{equation}
A = \frac{T}{2}
\end{equation}
$$

and suggests a value of $T=8$. This provides for the exact same value as Øye, but note that he uses slightly different notation than what is presented here.

Larsen used a dimensionalized and non-dimensionalized frequencies:

$$
\begin{equation}
A = \frac{1}{2\omega}
\end{equation}
$$

where 
$$
\begin{equation*}
\hat{\omega} = \frac{\omega c}{2v}
\end{equation*}
$$


#### Static Separation Point
Øye assumes a linear interpolation between the inviscid and static lifts, that interpolation percent he deems as the static separation point:

$$
\begin{equation}
f^s(\alpha) = \frac{C_l^{st} - C_l^{fs}}{C_l^{inv} - C_l^{fs}}
\end{equation}
$$


The static attatchment degree $f$ is used to describe the separation point of the fluid. The range of values for $f$ is in between $0$ and $1$: where $1$ would indicate fully attatched flow and $0$ represents fully separated flow. If $f$ were to equal a value of $0.75$, then the fluid would have a separation point that occurs $75\%$ the length of the chord away from the leading edge. 

For $f(\alpha)$, $C_l^{st}$ represents the  static lift coefficient; $C_l^{fs}$ is the lift for an airfoil that is under the condition of fully separated flow; and $C_l^{inv}$ is the lift for an airfoil in inviscid flow. It should also be mentioned that all of these coefficient of lift values are functions of angle of attack that range from $\alpha_{0} \leq \alpha \leq \alpha_{sep}$: where $\alpha_0$ is the angle of attack that produces zero lift and $\alpha_{sep}$ is the angle of attack where full separation occurs following static stall.

Because $C_l^{inv}$ is the lift under inviscid flow, it has a completely linear relationship with the angle of attack. The slope for the line of $C_l^{inv}$ is equal to the slope of $C_l^{st}$ at $\alpha_0$ (i.e. $\frac{\partial C_l^{st}}{\partial \alpha}\big|_{\alpha_0}$). This initial slope is approximately equal to $2\pi$, which allows $C_l^{inv}$ to be represented by the equation:

$$C_l^{inv} \approx 2\pi(\alpha - \alpha_0)$$

Note that here we default to $2\pi$, but ultimately the slope is set by the user. 

When Hansen describes this method in his textbook, he uses the same approach he used in his dynamic stall model:

$$
\begin{equation}
f^s_{H}\left(\alpha\right) = \left(2\sqrt{\frac{C_l^{st}}{C_l^{inv}}}\right)^2
\end{equation}
$$

##### Fully Separated Lift Coefficient
As was the case with finding the separation point, there are multiple ways to find the fully separated coefficient of lift. In Øye's original paper, he uses quadratic interpolation between $C^{s}_{l}(\alpha_0)$ and $C_l^{s}(\alpha_{sep})$, with a derivative constraint at $\alpha_0$. Thus the fully separated lift can be found by:

$$
\begin{equation}
 C_l^{fs}(\alpha) = a\alpha^2 + b\alpha + c 
\end{equation}
$$

where: 

$$
\begin{aligned}
a = & \frac{\frac{\partial C_l}{\partial \alpha}\big|_{\alpha_0}(\alpha_0-\alpha_{sep}) + C_{l_{sep}}}{(\alpha_0-\alpha_{sep})^2} \\
b = & \frac{-\alpha_0^2\frac{\partial C_l}{\partial \alpha}\big|_{\alpha_0} -2\alpha_0C_{l_{sep}} + \frac{\partial C_l}{\partial \alpha}\big|_{\alpha_0}\alpha_s^2}{(\alpha_0-\alpha_{sep})^2} \\
c = & \frac{\alpha_0(\alpha_0(C_{l_{sep}}+\frac{\partial C_l}{\partial \alpha}\big|_{\alpha_0}\alpha_s)-\frac{\partial C_l}{\partial \alpha}\big|_{\alpha_0}\alpha_s^2)}{(\alpha_0-\alpha_{sep})^2}
\end{aligned}
$$

Note that this approach is only valid between $\alpha_0$ and $\alpha_{sep}$. 

Since Hansen used a different approach to find the static separation point, he utilizes the linear interpolation between the inviscid and static lift coefficients to find the fully separated lift. 

$$
\begin{equation}
C_l^{fs}(\alpha) = \frac{C_l^{st} - f^s C_l^{inv}}{1-f^s}
\end{equation}
$$




Larsen and Faber use Hermite interpolation on the fully separated lift. Since the fully separated lift is unknown, they assume that the fully separated lift is equal to the static lift at $\alpha_0$ and $\alpha_{sep}$ and there derivatives are: 

 



$$
\frac{\partial C_l^{fs}}{\partial \alpha}|_{\alpha_0} = \frac{1}{2}(\frac{\partial C_l^{st}}{\partial \alpha}|_{\alpha_0})
$$

$$
\frac{\partial C_l^{fs}}{\partial \alpha}|_{\alpha_{sep}} = \frac{1}{12}(\frac{\partial C_l^{st}}{\partial \alpha}|_{\alpha_0})
$$

Using these derivatives, we find the Hermite interpolation to be:

$$
\begin{equation}
C_l^{fs}(\alpha) = t_0\left((\alpha_{sep}-\alpha_0)\frac{1}{2}\frac{\partial C_l^{st}}{\partial \alpha}\bigg|_{\alpha_0}\left(1 + t_0 \left( \frac{7}{6}t_1-1\right) \right)+C_{l_{sep}}^{st}t_0(1-2t_1) \right)
\end{equation}
$$

where:
$$
\begin{aligned}
t_0 = \frac{\alpha - \alpha_0}{\alpha_{sep}-\alpha_0} && t_1 = \frac{\alpha - \alpha_{sep}}{\alpha_{sep}-\alpha_0}
\end{aligned}
$$

As with Øye's original model, this model is only valid between $\alpha_0$ and $\alpha_{sep}$. For this approach it is traditionally assumed that $\alpha_{sep}=32^\circ$. Note that in our implementation of this function we assume that the fully separated lift converges to the static lift outside of this region.

Just as the separation point has an option of how to solve for it, DynamicStallModels.jl allows the user to toggle which method they want to use to find $C_l^{fs}$.




