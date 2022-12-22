# Beddoes-Leishman - Gonzalez
A translation of the OpenFAST v3.3.0 code. 

```math
\begin{align}
f''_{i-1} &= 0 \\
M &= \frac{U}{a} \\
\beta^2 &= 1 - M^2 \\
\beta &= \sqrt{1 - M^2} \\
C^{\text{circ}}_{n_\alpha} &= \frac{C_{n,\alpha}}{\beta} = \frac{1}{\beta}\frac{\partial C_n}{\partial \alpha} \\ 
\end{align}
```
- I think the second half of equation 5 is correct, but I should double check 
- If `FLookup`, then set $\eta_e=1.0$.
- There is some filtering going on. They have some values such as $\alpha_{filt}$ which they denote as `alpha_filt_curr` and `alpha_filt_minus1` and I'm not entirely sold on the notation. It's easy to mistake what you're working with. I'm renaming $\alpha_{filt}$ to $\theta$. In the future, I'd like to give the input $\alpha$ a special name, then just call everything else something make sense. 
- On the `FirstPass`, set $\alpha_{i-1}$ and $\theta_{i-1}$ to $\alpha_i$


```math
\begin{align}
\Delta s &= \frac{2U \Delta t}{c} \\
z &= 20 \\ 
\xi &= \frac{z}{\pi c}\max\left(1, U\right) \\
c_{lp} &= \exp\left(-2\pi\Delta t \xi \right) \\
\theta_i &= c_{lp}\theta_{i-1} + (1-c_lp)\alpha\\
\Delta \alpha_0 &= \theta_i - \alpha_0 \\
K_{\alpha_i} &= \frac{\theta_i - \theta_{i-1}}{\Delta t}
\end{align}
```

- On `FirstPass` set $K_{\alpha_{i-1}}= 0$. 


```math
\begin{align}
q_i = \frac{K_{\alpha_i}c}{U}
\end{align}
```

- On `FirstPass` set $q_{i-1}$ and $Q_{i-1}$ to $q_i$
- I did a name change on `q_f_curr` to $Q_i$

```math
\begin{align}
Q_i &= C_{lp}Q_{i-1} + (1-c_{lp})q_i \\
K_{\alpha, f} &= \frac{Q_i U}{c} \\
K_q &= \frac{q_i - q_{i-1}}{\Delta t} 
\end{align}
```
- Note that I don't know if $K_{\alpha, f}$ is a separate quantity, or a filtered quantity. 
- If `FirstPass` set $K_{q,f_{i-1}}=0$.

```math
\begin{align}
K_{q,f} &= c_{lp}K_{q,f_{i-1}} + (1-c_{lp})K_q \\
k_\alpha &= \frac{1}{(1 -M) + \frac{1}{2}\frac{\partial C_n}{\partial \alpha} M^2 \beta (A_1b_1 +A_2b_2)} \\
k_q &= \frac{1}{(1 -M) + \frac{\partial C_n}{\partial \alpha} M^2 \beta (A_1b_1 +A_2b_2)} \\
T_I &= \frac{c}{a} \\
T_\alpha &= \frac{3}{4}T_I k_\alpha \\
T_q &= \frac{3}{4}k_q \\
T_f &= \frac{T_{f0}}{\sigma_1} \\
T_{fc} &= \frac{T_{f0}}{\sigma_{1c}} \\
K'_{\alpha_i} &= f\left(\Delta t, T_\alpha, K'_{\alpha_{i-1}} \right)
\end{align}
```