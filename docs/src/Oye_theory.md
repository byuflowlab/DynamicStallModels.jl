#  Øye Theory

Øye's model is a linear interpolation of fully attached lift (inviscid), and fully detached lift based on a dynamic degree of attachment. The dynamic attachment degree is a simple first-order filter on the separation point on the suction side of an airfoil. This dynamic degree of attachment is the Øye model's only state, which makes it the simplest dynamic stall model to implement.

The state is:

| Variable | Name | Comments |
| -------- | -------------------------- | ----- |
| $\dot{f}^d$ | Dynamic Attachment Degree | This differential equation is also utilized in the Riso model and Beddoes-Leishman model |

## State Rate Equation

$$
\begin{equation}
\dot{f}^d = -\tau (f^d- f(\alpha))
\end{equation}
$$

As stated, $f^d$ is the only state rate equation in this model. $\tau$ is the time constant, and $f^s(\alpha)$ is the position of the fluid's separation point as a function of the angle of attack.

Additionally, it should be noted that variables with a dot above them, such as $\dot{f}^d$, are derivatives with respect to time.

The above equation can be evaluated using a prefered ODE solver. Once the solution is discovered, it can be implemented into the following equation to solve for the dynamic lift.

$$
\begin{equation}
C_{L}(t) = f^{d}(t)C_l^{inv} + (1-f^{d}(t))C_l^{fs}
\end{equation}
$$

The indicial approach simply solves equation 1 by assuming that the inflow velocity and angle of attack do not vary within the time step, and they use an integrating factor:

$$
\begin{equation}
f^{d}_{i+1} = f^{s}(\alpha_i) + (f^{d}_i - f^{s}(\alpha_i))\hspace{3pt}e^{-\frac{dt}{\tau}}
\end{equation}
$$


## Explanation of Variables
#### Time Constant

$$
\begin{equation}
\tau = \frac{V}{A c}
\end{equation}
$$

The time constant is experimentally determined, and provides the needed information for how long it takes to come to steady state. Each author tends to implement the time constant in a different way. Here we rely on Hansen's approach (as seen above). 

We need to show each of the author's equations for the time constant and how to switch between (we don't have to show our work though.). 
Faber's Approach
The constant $\tau$ is said to be around $\tau \approx 8$.

Larsen's Approach


#### Static Separation Point

$$
\begin{equation}
f^s(\alpha) = \frac{C_l^{s} - C_l^{fs}}{C_l^{inv} - C_l^{fs}}
\end{equation}
$$

There are a few common ways to calculate the separation point. The way shown above is how Øye, Larsen, and Faber obtain this value. However, in Hansen's version of the Øye method, he uses a noticeably different approach:

$$
\begin{equation}
f^s_{H}\left(\alpha\right) = \left(2\sqrt{\frac{C_l^{st}}{C_l^{inv}}}\right)^2
\end{equation}
$$

Note that we notate the separation point from Hansen's paper with the subscript $H$, so the two versions of the separation point can be distinguished.

In DynamicStallModels.jl, either one of these ways of finding the separation point can be implemented with a simple toggle. 

The static attatchment degree $f$ is used to describe the separation point of the fluid. The range of values for $f$ is in between $0$ and $1$: where $1$ would indicate fully attatched flow and $0$ represents fully separated flow. If $f$ were to equal a value of $0.75$, then the fluid would have a separation point that occurs $75\%$ the length of the chord away from the leading edge. 

For $f(\alpha)$, $C_l^{s}$ represents the  static lift coefficient; $C_l^{fs}$ is the lift for an airfoil that is under the condition of fully separated flow; and $C_l^{inv}$ is the lift for an airfoil in inviscid flow. It should also be mentioned that all of these coefficient of lift values are functions of angle of attack that range from $\alpha_{0} \leq \alpha \leq \alpha_{sep}$: where $\alpha_0$ is the angle of attack that produces zero lift and $\alpha_{sep}$ is the angle of attack where full separation occurs following static stall.

Because $C_l^{inv}$ is the lift under inviscid flow, it has a completely linear relationship with the angle of attack. The slope for the line of $C_l^{inv}$ is equal to the slope of $C_l^{s}$ at $\alpha_0$ (i.e. $\frac{\partial C_l^{s}}{\partial \alpha}\big|_{\alpha_0}$). This initial slope is approximately equal to $2\pi$, which allows $C_l^{i}$ to be represented by the equation:

$$C_l^{inv} \approx 2\pi(\alpha - \alpha_0)$$

Note that here we default to $2\pi$, but ultimately the slope is set by the user. 

##### Fully Separated Lift Coefficient
As was the case with finding the separation point, there are multiple ways to find the fully separated coefficient of lift. In Øye's original paper, he uses quadratic interpolation between $C^{s}_{l}(\alpha_0)$ and $C_l^{s}(\alpha_{sep})$. Additionally it has a derivative constraint at $\alpha_0$. 

Hansen finds this coefficient differently:

$$
\begin{equation}
C_l^{fs} = \frac{C_l^{st} - f^s C_l^{inv}}{1-f^s}
\end{equation}
$$




However, both Larsen and Faber utilizes not only the derivative at $\alpha_0$ but also at $\alpha_{sep}$. Faber notes in his paper that an angle of attack of $32 \degree$ is a good assumption of where full separation occurs. They then use Hermite Interpolation to create a polynomial that gives the fully separated lift. Including the extra derivative at $\alpha_{sep}$ in the curve fitting process adds extra robustness to this polar. This is shown here:

- [] I want to combine all of the Hermite interpolation stuff into a single spot, we don't need to go so far into detail. 

$$
\frac{\partial C_l^{fs}}{\partial \alpha}|_{\alpha_0} = \frac{1}{2}(\frac{\partial C_l^{st}}{\partial \alpha}|_{\alpha_0})
$$

$$
\frac{\partial C_l^{fs}}{\partial \alpha}|_{\alpha_{sep}} = \frac{1}{12}(\frac{\partial C_l^{st}}{\partial \alpha}|_{\alpha_0})
$$


Just as the separation point has an option of how to solve for it, DynamicStallModels.jl allows the user to toggle which method they want to use to find $C_l^{fs}$.

#### Hermite Interpolation

Hermite Interpolation is used to find the intermediary values of $C_l^{fs}$ between $\alpha_0$ and $\alpha_{sep}$. Faber's paper gives a good description of this process, and a similar format to that paper will be used to describe Hermite Interpolation.

For this interpolation method, the starting point $(x_0,y_0)$ and the ending point $(x_1,y_1)$ must be known. Additionally, $\frac{dy}{dx}|_{0}$ and $\frac{dy}{dx}|_{1}$ must be given. With these values, the interpolation process is as follows:

$$t_0 = \frac{x-x_0}{x_1-x_0} ~~ ~~ \mathrm{and} ~~ ~~ t_1 = \frac{x-x_1}{x_1-x_0}$$

$$y(x) = \frac{dy}{dx}|_{0} ~(x_1-x_0)(t_0 + (t_1-1)t_0^2)~ + ~ \frac{dy}{dx}|_{1} ~ (x_1-x_0)t_0^2t_1 ~ + ~ y_1t_0^2(1-2t_1)$$

The values above can be switched to match what is needed for the Øye method:

$$
\begin{align}
1
\end{align}
% y(x) = C_l^{fs} ~,~x = \alpha ~,~ x_0 = \alpha_0 ~,~ x_1 = \alpha_{sep} ~,~ y_0 = 0 ~,~ y_1 = C_l^{st}(\alpha_{sep}) ~,~ \frac{dy}{dx}|_{0}=\frac{1}{2}(\frac{\partial C_l^{st}}{\partial \alpha}|_{\alpha_0}) ~,~ \frac{dy}{dx}|_{1} = \frac{1}{12}(\frac{\partial C_l^{st}}{\partial \alpha}|_{\alpha_0})
$$

The fully separated lift coefficients can now be determined at angle of attack values between $\alpha_0$ and $\alpha_{sep}$.




