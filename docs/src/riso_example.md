# Risø

This reproduces the former introductory pitching-flat-plate inputs using the
repaired **continuous aggregate wrapper**. Use the separate example environment
from [Getting Started](@ref), starting at the repository root.

The wrapper requires `A`, `b` and `T` explicitly. Its `T=[Tp,Tf]` entries are
**seconds**, matching the existing rate equations. Do not substitute
nondimensional time constants from a paper without the appropriate conversion.
These historical example values are retained for equation reproduction; this
page does not establish their empirical calibration or separated-flow validity.

```julia
using DynamicStallModels, OrdinaryDiffEqTsit5
ds = DynamicStallModels
u = c = 1.0
omega = 0.2 * 2u/c
U(t) = u
Udot(t) = 0.0
alpha(t) = (2sin(omega*t) + 5)*(pi/180)
alphadot(t) = 2omega*cos(omega*t)*(pi/180)
a = collect(-pi:0.01:pi)
polar = hcat(a, 2pi .* a, zeros(length(a)))
af = make_airfoil(polar, ds.NoModel(); sfun=ds.RSP())
A = [0.165, 0.335]
b = [0.0455, 0.3000]
T = [3.0, 6.0]
model = Riso([af]; A, b, T)
x0 = [0.0, 0.0, 1.0, 0.0]
parameters = [U, Udot, alpha, alphadot, c]
problem = ds.SciMLBase.ODEProblem{true}(model, x0, (0.0, 80.0), parameters)
times = collect(0.0:0.1:80.0)
solution = solve(problem, Tsit5(); dtmax=0.1, saveat=times)
@assert ds.SciMLBase.successful_retcode(solution)
Cl, Cd, t = parsesolution(model, solution, parameters)
@assert size(Cl) == size(Cd) == (length(t), 1)
@assert all(isfinite, Cl) && all(isfinite, Cd)
```

The four state blocks are kept separate for every station. For different
response constants use separate aggregate models; the constants above are shared
by the stations within this wrapper. The underlying scalar rate functions and
historical indicial recurrence are separate interfaces. Passing `Discrete()` to
this aggregate wrapper is rejected instead of returning rates as if they were
discrete updated states.

For this linear polar, attachment forcing is one. The two wake states have
closed-form sinusoidal responses, giving an independent lift check:

```julia
Tu = c/(2u)
function wake_exact(j, time)
    rate = b[j]/Tu
    particular(s) = A[j]*5pi/180 + A[j]*2pi/180 *
        (rate^2*sin(omega*s)-rate*omega*cos(omega*s))/(rate^2+omega^2)
    particular(time)-particular(0.0)*exp(-rate*time)
end
ae = alpha.(t).*(1-sum(A)) .+ wake_exact.(1,t) .+ wake_exact.(2,t)
attachment = 1 .- exp.(-t/T[2])
exact_Cl = 0.5af.dcldalpha .* (ae.-af.alpha0).*(1 .+ attachment) .+ pi*Tu.*alphadot.(t)
max_Cl_error = maximum(abs.(Cl[:,1].-exact_Cl))
@show max_Cl_error
@assert max_Cl_error < 1e-6

using Plots
fig = plot(t, Cl[:,1]; label="Risø", xlabel="Time (s)", ylabel="Cl")
plot!(fig, t, exact_Cl; label="Analytical attached-flow solution", linestyle=:dash)
plot_path = joinpath(mktempdir(), "riso_attached.png")
savefig(fig, plot_path)
@assert isfile(plot_path)
```

For coefficient extraction without an aggregate model, supply the response
amplitudes: `ds.riso_coefficients(x,U,alpha,alphadot,c,af; A=A)`.
Residuals similarly require explicit constants:
`ds.riso_residual(dx,x,[U,Udot,alpha,alphadot],c,t,af; A,b,T)`.
The legacy `static_riso_dae!` takes the Airfoil as its sixth argument and returns
`dx-rates`, preserving the supplied derivative vector. Use zero `dx` for an
equilibrium residual. None of these repairs establish experimental stall agreement.
