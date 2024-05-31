# Recursive Solution to Duhamel's Integral

First we start with Duhamel's integral:

$$
\begin{equation}
    C_l(t) = \frac{\partial C_l}{\partial \alpha} \left[\alpha(t_0)\phi(s) + \int\limits_{s_0}^s \frac{d\alpha}{dt}(\sigma) \phi(s - \sigma) d\sigma \right] = \frac{\partial C_l}{\partial \alpha} \alpha_e(t)
\end{equation}
$$

Where $\alpha_e$ is the equivalent angle of attack. 

Assume there exists a solution such that

$$
\begin{equation}
    \phi(s) = 1 - A_1e^{-b_1 s} - A_2e^{-b_2 s}
\end{equation}
$$

Then the equivalent angle of attack becomes: 

$$
\begin{aligned}
\alpha_e(s) = \alpha(s_0)\left(1 - A_1e^{-b_1 s} - A_2e^{-b_2 s} \right)\\ + \int\limits_{s_0}^s \frac{d\alpha}{dt}(\sigma) \left(1 - A_1e^{-b_1 (s-\sigma)} - A_2e^{-b_2 (s-\sigma)} \right) d\sigma \\
= \alpha(s_0) - A_1 \alpha(s_0)e^{-b_1 s} - A_2 \alpha(s_0)e^{-b_2 s}\\ + \int\limits_{s_0}^s d\alpha(s) -  A_1 \int\limits_{s_0}^s\frac{d \alpha}{d s}(\sigma)e^{-b_1 (s-\sigma)} d\sigma \\ -  A_2 \int\limits_{s_0}^s\frac{d \alpha}{d s}(\sigma)e^{-b_2 (s-\sigma)} d\sigma
\end{aligned}
$$

Note that the terms $A_1 \alpha(s_0) e^{-b_1 s}$ and $A_2 \alpha(s_0) e^{-b_2 s}$ will quickly drop to zero (assume them zero). Then the equivalent angle of attack becomes:

$$
\begin{equation}
\alpha_e(s) = \alpha(s) - X(s) - Y(s)
\end{equation}
$$

where 

$$
\begin{equation}
X(s) = A_1\int\limits_{s_0}^s\frac{d \alpha}{d s}(\sigma)e^{-b_1 (s-\sigma)} d\sigma
\end{equation}
$$

and

$$
\begin{equation}
Y(s) = A_2 \int\limits_{s_0}^s\frac{d \alpha}{d s}(\sigma)e^{-b_2 (s-\sigma)} d\sigma
\end{equation}
$$

If we consider just $X(s)$, then the same goes for $Y(s)$. So first we assume that $s_0 = 0$, then let's assume we can find the value after some time $s$ (So basically what we have above, just replace $s_0$ with 0.). Then if we want to find it some time ($s + \Delta s$) after we have

$$
\begin{equation}
    X(s + \Delta s) = A_1\int\limits_{0}^{s+\Delta s}\frac{d \alpha}{d s}(\sigma)e^{-b_1 (s+\Delta s-\sigma)} d\sigma \\
\end{equation}
$$

We can split this integral up into: 

$$
\begin{aligned}
X(s + \Delta s) =A_1 e^{-b_1 \Delta s}\int\limits_{0}^s\frac{d \alpha}{d s}(\sigma)e^{-b_1 (s-\sigma)} d\sigma \\ + A_1\int\limits_{s}^{s+\Delta s}\frac{d \alpha}{d s}(\sigma)e^{-b_1 (s+\Delta s-\sigma)} d\sigma
\end{aligned}
$$

But we see the original $X(s)$ is that first integrand, so we have: 

$$
\begin{equation}
X(s + \Delta s) =X(s) e^{-b_1 \Delta s} + A_1\int\limits_{s}^{s+\Delta s}\frac{d \alpha}{d s}(\sigma)e^{-b_1 (s+\Delta s-\sigma)} d\sigma
\end{equation}
$$

Which is really just the previous state (the first term) decayed across time, with an update (the second term). We'll re-label the second term as $I$. If we consider use the following second-order backwards difference for $\frac{d \alpha}{ds}$, then we will get an integral that we can evaluate the update term exactly:

$$
\begin{equation}
    \frac{d \alpha}{ds}\bigg|_{s+\Delta s} = \frac{3 \alpha(s + \Delta s) - 4\alpha (s) + \alpha (s - \Delta s)}{2\Delta s}
\end{equation}
$$

then that update term is evaluated as follows: 

$$
\begin{aligned}
I = A_1\int\limits_{s}^{s+\Delta s}\frac{d \alpha}{d s}(\sigma)e^{-b_1 (s+\Delta s-\sigma)} d\sigma \\
= A_1e^{-b_1 (s+\Delta s)}\int\limits_{s}^{s+\Delta s}\frac{3 \alpha(s + \Delta s) - 4\alpha (s) + \alpha (s - \Delta s)}{2\Delta s}e^{b_1\sigma} d\sigma \\
A_1 \left( \frac{\Delta \alpha (s + \Delta s)}{\Delta s}\right)\left(\frac{1 - e^{-b_1 \Delta s}}{b_1} \right)
\end{aligned}
$$

If the product $b_1 \Delta s$ is small, we can neglect higher order terms, so 

$$
\begin{equation}
\frac{1 - e^{-b_1 \Delta s}}{b_1} \approx \Delta s.
\end{equation}
$$

Then the update term becomes:

$$
\begin{equation}
I = A_1 \left( \frac{\Delta \alpha (s + \Delta s)}{\Delta s}\right) \Delta s = A_1 \Delta \alpha (s + \Delta s)
\end{equation}
$$

This turns out to be equivalent to assuming $\frac{d \alpha}{d s}$ is constant across the time step and using the left-hand rectangular rule for evaluating the update integral. :|

But this results in a simple update formula:

$$
\begin{equation}
X(s+ \Delta s) = X(s) e^{-b_1 \Delta s} +  A_1 \Delta \alpha (s + \Delta s)
\end{equation}
$$

Generally to obtain errors of less than 5%, each of the products $b_i \Delta s$ must be less than 0.05. 

If you use the midpoint rule instead of the left-hand rule, then the increment is evaluated as follows: 

$$
\begin{equation}
I = A_1 \int\limits_s^{s+\Delta s} \frac{d \alpha}{ds} e^{-b_1(s + \Delta s - \sigma)}d\sigma 
\end{equation}
$$

and if you assumed that the step change of $\alpha$ is constant across the time step, then 
$$
\begin{equation}
I = A_1 \left[\frac{\Delta \alpha}{\Delta s} e^{-b_1(s + \Delta s - (s + \Delta s/2))} \Delta s \right] = A_1 \Delta \alpha e^{-b_1 \Delta s/2}
\end{equation}
$$

Which means that for the total time step update you have: 

$$
\begin{equation}
X(s+ \Delta s) = X(s) e^{-b_1 \Delta s} +  A_1 \Delta \alpha e^{-b_1 \Delta s/2}
\end{equation}
$$

Generally to obtain errors of less than 1%, each of the products $b_i \Delta s$ must be less than 0.25. 

The midpoint rule is what is used in the Beddoes-Leishman model. 
