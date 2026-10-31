# Øye Theory

The supported lift-only Øye model blends inviscid and fully separated lift
through one delayed attachment state f, normally between zero and one:

```math
\frac{df}{dt}=\frac{f^s(\alpha)-f}{\tau},\qquad
\tau=\frac{A c}{U},\qquad
C_l=f C_l^{inv}+(1-f)C_l^{fs}.
```

Here c is chord in metres, U is relative-speed magnitude in m/s, A is a
positive empirical coefficient (commonly 4), and tau is a **time in seconds**.
The inverse time is U/(A*c); it must not also be called the dimensional time
constant. At U=0 this implementation holds the state fixed. Negative relative
speed and nonpositive chord/A are rejected.

For held inputs over a time step, the exact scalar update is

```math
f_{i+1}=f^s(\alpha_{i+1})+
       \left[f_i-f^s(\alpha_{i+1})\right]\exp(-\Delta t/\tau_{i+1}).
```

The discrete solver uses right-end samples; the continuous solver evaluates
input functions at the current integration time. These approximations converge
with sampling refinement for smooth forcing, but are not identical at finite dt.

## Static and fully separated lift

The inviscid coefficient uses the supplied/fitted slope and zero-lift angle:

```math
C_l^{inv}=C_{l\alpha}(\alpha-\alpha_0).
```

For the Hansen (`version=1`, `RSP()`) convention, away from zero inviscid lift,

```math
f^s=\left(2\sqrt{|C_l^{st}/C_l^{inv}|}-1\right)^2,\qquad
C_l^{fs}=\frac{C_l^{st}-f^s C_l^{inv}}{1-f^s}.
```

Attachment is bounded to [0,1] and the RSP function applies its static-angle
limits. At f^s=1 the separated-lift quotient is 0/0 for an exactly attached
polar. Along the Hansen relation its limit is Cinv/2; the implementation uses
that finite convention. For an arbitrary user separation function returning 1,
Cfs is unidentifiable from equilibrium data and the same convention is used.

The Faber/Larsen (`version=2`, `LSP()`) option retains the positive-incidence
Hermite polynomial. It joins the static coefficient at alpha0 and positive
alpha_sep with endpoint slopes C_lalpha/2 and C_lalpha/12. Outside that interval
it returns the static coefficient; no new negative-stall fit is supplied.
The existing 0.35 near-equality rule in LSP is a heuristic, not a universal
physical constant or proof of static-polar recovery everywhere.

This model represents attachment lag, not explicit leading-edge-vortex
advection/shedding. Only Cl is modelled here; Cd and Cm output columns are zero.
Refer to [Øye](@ref) for runnable examples and the [Audit and Fix Plan](@ref)
for the scope of verification versus scientific validation.
