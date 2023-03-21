#  Øye Theory

Øye's model is a linear interpolation of fully attached lift (inviscid), and fully detached lift based on a dynamic degree of attachment. The dynamic attachment degree is a simple first-order filter on the separation point on the suction side of an airfoil. This dynamic degree of attachment is the Øye model's only state, which makes it the simplest dynamic stall model to implement.

The state is:

| Variable | Name | Comments |
| -------- | -------------------------- | ----- |
| $\dot{f_{dyn}}$ | Dynamic Attatchment Degree | This differential equation is also utilized in the Riso model and Beddoes-Leishman model |

## State Rate Equation

 $$\dot{f_{dyn}} = -T_f (f_{dyn}- f(\alpha))$$

As stated, $f_{dyn}$ is the only state rate equation in this model. $T_f$ is called the time scale parameter, and $f(\alpha)$ is the position of the fluid seperation point as a function of the angle of attack.

Additionally, it should be noted that variables with a dot above them, such as $\dot{f_{dyn}}$, are derivatives with respect to time.

## Time Scale Parameter

$$T_f = \frac{2V}{\tau c}$$

$T_f$ represents the time scale parameter for dynamic stall. In a way, it describes the combined rate for the development of the vortex that forms during dynamic stall and its shedding off of the airfoil. Also, $V$ is the incoming flow velocity, and $c$ is the chord length of the airfoil. The constant $\tau$ is said to be around $\tau \approx 8$.
## Separation Point as a Function of Alpha

$$f(\alpha) = \frac{C_l^{st} - C_l^{fs}}{C_l^{inv} - C_l^{fs}}$$

There are a few accepted ways to calculate the seperation point. The way shown above is how Øye, Larsen, and Faber obtain this value. However, in Hansen's version of the Øye method, he uses a noticeably different approach:

$$f_{H} = (2\sqrt{\frac{C_l^{st}}{C_l^{inv}}})^2$$

Note that the seperation point in Hansen's paper is denoted with the subscript $H$, so the two versions of the seperation point can be distinguished.

In DynamicStallModels.jl, either one of these ways of finding the seperation point can be implemented with a simple toggle. 

The variable $f$, the static attatchment degree, is used to describe the seperation point of the fluid. The range of values for $f$ is in between $0$ and $1$: where $1$ would indicate fully attatched flow and $0$ represents fully seperated flow. If $f$ were to equal a value of $0.75$, then the fluid would have a seperation point that occurs $75\%$ the length of the chord away from the leading edge. 

For $f(\alpha)$, $C_l^{st}$ represents the  static lift coefficient; $C_l^{fs}$ is the lift for an airfoil that is under the condition of fully seperated flow; and $C_l^{inv}$ is the lift for an airfoil in inviscid flow. It should also be mentioned that all of these coefficient of lift values are functions of angle of attack that range from $\alpha_{0} \leq \alpha \leq \alpha_{sep}$: where $\alpha_0$ is the angle of attack that produces zero lift and $\alpha_{sep}$ is the angle of attack where full seperation occurs following static stall.

Because $C_l^{inv}$ is the lift under inviscid flow, it has a completely linear relationship with the angle of attack. The slope for the line of $C_l^{inv}$ is equal to the slope of $C_l^{st}$ at $\alpha_0$ (i.e. $\frac{\partial C_l^{st}}{\partial \alpha}|_{\alpha_0}$). This initial slope is approximately equal to $2\pi$, which allows $C_l^{inv}$ to be represented by the equation:

$$C_l^{inv} = 2\pi(\alpha - \alpha_0)$$

As was the case with finding the seperation point, there are multiple ways to find the fully separated coefficient of lift.

In Øye's original paper, he simply matches a parabola that has the derivative constraint at $\alpha_0$ and passes through the point $(\alpha_0, 0)$ and $(\alpha_{sep}, C_l^{st}(\alpha_{sep}))$. However, both Larsen and Faber utilizes not only the derivative at $\alpha_0$ but also at $\alpha_{sep}$. Faber notes in his paper that an angle of attack of $32 \degree$ is a good assumption of where full seperation occurs. They then use Hermite Interpolation to create a polynomial that gives the fully separated lift. Including the extra derivative at $\alpha_{sep}$ in the curve fitting process adds extra robustness to this polar.

The accepted values of these derivatives from Øye, Larsen, and Faber are given below:

$$\frac{\partial C_l^{fs}}{\partial \alpha}|_{\alpha_0} = \frac{1}{2}(\frac{\partial C_l^{st}}{\partial \alpha}|_{\alpha_0})$$
$$\frac{\partial C_l^{fs}}{\partial \alpha}|_{\alpha_{sep}} = \frac{1}{12}(\frac{\partial C_l^{st}}{\partial \alpha}|_{\alpha_0})$$


Hansen finds this coefficient differently:

$$C_l^{fs} = \frac{C_l^{st} - f_H C_l^{inv}}{1-f_H}$$

Just as the seperation point has an option of how to solve for it, DynamicStallModels.jl allows the user to toggle which method they want to use to find $C_l^{fs}$.

## Hermite Interpolation

Hermite Interpolation is used to find the intermediary values of $C_l^{fs}$ between $\alpha_0$ and $\alpha_{sep}$. Faber's paper gives a good description of this process, and a similar format to that paper will be used to describe Hermite Interpolation.

For this interpolation method, the starting point $(x_0,y_0)$ and the ending point $(x_1,y_1)$ must be known. Additionally, $\frac{dy}{dx}|_{0}$ and $\frac{dy}{dx}|_{1}$ must be given. With these values, the interpolation process is as follows:

$$t_0 = \frac{x-x_0}{x_1-x_0} ~~ ~~ \mathrm{and} ~~ ~~ t_1 = \frac{x-x_1}{x_1-x_0}$$

$$y(x) = \frac{dy}{dx}|_{0} ~(x_1-x_0)(t_0 + (t_1-1)t_0^2)~ + ~ \frac{dy}{dx}|_{1} ~ (x_1-x_0)t_0^2t_1 ~ + ~ y_1t_0^2(1-2t_1)$$

The values above can be switched to match what is needed for the Øye method:

$$y(x) = C_l^{fs} ~,~x = \alpha ~,~ x_0 = \alpha_0 ~,~ x_1 = \alpha_{sep} ~,~ y_0 = 0 ~,~ y_1 = C_l^{st}(\alpha_{sep}) ~,~ \frac{dy}{dx}|_{0}=\frac{1}{2}(\frac{\partial C_l^{st}}{\partial \alpha}|_{\alpha_0}) ~,~ \frac{dy}{dx}|_{1} = \frac{1}{12}(\frac{\partial C_l^{st}}{\partial \alpha}|_{\alpha_0})$$

The fully seperated lift coefficients can now be determined at angle of attack values between $\alpha_0$ and $\alpha_{sep}$.

## Solving The State Rate Equation

$$\dot{f_{dyn}} = -T_f (f_{dyn} - f(\alpha))$$

Now that values for $T_f$ and $f(\alpha)$ can be solved for, this equation can be evaluated using a prefered ODE solver. 


$$C_L(t) = f_{dyn}(t)C_l^{inv} + (1-f_{dyn}(t))C_l^{fs}$$


