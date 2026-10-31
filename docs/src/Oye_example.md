# Øye

Øye delays one attachment state f. These examples use the supported **lift-only**
mode: `cflag=1`. Outputs are Cl, zero Cd and zero Cm; zero columns are not drag
or moment predictions. `cflag=2` now raises a clear error because its previous
normal-force implementation was empty. No new force-closure assumption is made.

Use the example environment from [Getting Started](@ref) and start from the
repository root. Copy the blocks in order. The stored Faber polar's first column
already spans -pi to pi in radians; do not convert it again.

## Discrete solve

```julia
using DynamicStallModels, DelimitedFiles
ds = DynamicStallModels
polar = readdlm(joinpath(pwd(), "data", "NACA_0015_Faber.csv"), ',', Float64)
model = Oye(Discrete(), 1, 2, 4.0)
af = ds.make_airfoil(polar, model; sfun=ds.LSP())
c = 0.55
speed = 0.11 * 343.0
omega = 0.051 * 2speed / c
alpha(t) = deg2rad(10 + 10sin(omega*t))
t = collect(range(0.0, 2.0; length=2001))
U = fill(speed, length(t))
states, loads = solve_discrete([af], t, U, alpha.(t); cvec=[c])
@assert size(states) == (length(t), 1) && size(loads) == (length(t), 3)
@assert all(isfinite, loads) && all(f -> 0 <= f <= 1, states)
@assert all(iszero, loads[:, 2:3])
@assert maximum(abs, loads[:, 1]) < 3.0 # broad plausibility check, not validation
```

Version 2 uses the existing Faber/Larsen positive-incidence Hermite construction
and its LSP separation function. Outside alpha0–positive separation angle the
fully separated lift uses the static polar; this is not a calibrated symmetric
negative-stall extension. The existing LSP near-equality heuristic is retained.

`solve_discrete` initializes f from the first static separation value and uses
the right-end input for each held step. The dimensional response time is
tau=A*c/U. `solve_indicial` is retained as a compatibility alias.

## Continuous solve

The same physical model can be integrated with an ODE solver. Geometry is
separate from the four environment entries; omitting `cvec` would use 1 m.
The derivatives are unused by Øye, so they may be numbers rather than functions.

```julia
using OrdinaryDiffEqTsit5
continuous = Oye(Continuous(), 1, 2, 4.0)
afc = ds.make_airfoil(copy(polar), continuous; sfun=ds.LSP())
environment = [t -> speed, 0.0, alpha, 0.0]
x0 = [states[1, 1]]
problem = ODEProblem([afc], x0, (first(t), last(t)), environment; cvec=[c])
solution = solve(problem, Tsit5(); reltol=1e-9, abstol=1e-11,
                 saveat=t, dtmax=minimum(diff(t)))
@assert ds.SciMLBase.successful_retcode(solution)
@assert all(isfinite, Array(solution))
answer = parsesolution(continuous, afc, solution, environment, [c, 0.2])
@assert size(answer) == (2, length(t))
```

`answer` has angle in radians in row 1 and Cl in row 2. It uses the same load
routine as the discrete solve. Under time-varying forcing the held-input
discrete solution need not equal the continuously forced ODE solution at finite
dt; refine the sampling before interpreting the difference as model error.
`saveat` only selects output times: it does not limit internal integration steps.
The `dtmax` limit prevents an adaptive solver from skipping short changes in the
polar-based separation forcing during a nearly constant part of the response.

```julia
using Plots
fig = plot(rad2deg.(alpha.(t)), loads[:, 1]; label="Discrete", xlabel="Angle (deg)", ylabel="Cl")
plot!(fig, rad2deg.(answer[1, :]), answer[2, :]; label="Continuous", linestyle=:dash)
plot_path = joinpath(mktempdir(), "oye_lift.png")
savefig(fig, plot_path)
@assert isfile(plot_path)
@show plot_path
```

For the Hansen fully separated lift, use `Oye(Discrete(),1,1,4.0)` with
`sfun=ds.RSP()`. Its fully attached removable singularity uses the Hansen limit
Cinv/2. See [Øye Theory](@ref). Numerical verification tests are active in
`test/Oye_Test.jl`; they do not establish experimental agreement or full dynamic
stall validation of these formulations.
