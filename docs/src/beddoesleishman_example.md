# Beddoes–Leishman

## AeroDyn/Gonzalez discrete implementation

Begin with [Getting Started](@ref): it is a complete BL-ADG run, including
environment setup, units, data loading, finite-output checks and labelled
comparison plots. Its code is executed by `test/run_documentation_examples.jl`.

The original and continuous state-space implementations remain under
construction. The usable path here is `Discrete()` version 3 (Gonzalez),
created by the pinned OpenFASTTools helper. This is OpenFAST UA_Mod=2, not 3.

The following compact reproduction is independent of the Getting Started
session. Run it from the repository root in `--project=examples`:

```julia
using DynamicStallModels, OpenFASTTools, DelimitedFiles
repo = pwd()
input = OpenFASTTools.read_airfoilinput(joinpath(repo, "data", "DU21_A17.dat"))
af, xcp = OpenFASTTools.make_dsairfoil(input)
raw = readdlm(joinpath(repo, "data", "aerodynout_fordynamicstall.csv"), ',')
names = String.(raw[1, :]); names[1] = "Time"
values = Float64.(raw[3:end, :])
outs = Dict(names[i] => values[:, i] for i in eachindex(names))
t = outs["Time"]
U = hypot.(outs["AB1N011Vx"], outs["AB1N011Vy"])
alpha = deg2rad.(outs["AB1N011Alpha"])
states, loads = solve_discrete([af], t, U, alpha; cvec=[1.0], xcpvec=[xcp])
@assert size(states) == (100, 32) && size(loads) == (100, 3)
@assert all(isfinite, states) && all(isfinite, loads)
Cl, Cd, Cm = loads[:, 1], loads[:, 2], loads[:, 3]
Cn = Cl .* cos.(alpha) .+ Cd .* sin.(alpha)
Ct = Cl .* sin.(alpha) .- Cd .* cos.(alpha)
```

Each station contributes **Cl, Cd, Cm**, not five output columns. For several
stations the columns are `[Cl1,Cd1,Cm1,Cl2,Cd2,Cm2,...]`. The input histories in
this solver are shared across the stations; chord/xcp are per station. The first
row is a static initialization; later rows use dynamic coefficients.

State columns relevant to investigation are 1=unfiltered alpha (radians),
23=delayed normal separation fpp, 26=vortex driving source Cv, 27=vortex normal
state Nv, 28=convective vortex clock tauv, and 29–32=flags/first-pass state.
A nonzero Nv is not automatically a nonzero applied vortex force: clock and
cutout blending also matter.

## What the reference establishes

The supplied low-incidence CSV is rounded and its original OpenFAST binary,
UA variant and full initialization were not recovered. Historical consumers
used AeroDyn-original; the later Gonzalez example reuses that file. Current
modified-NREL5MW blade/driver parameters come from a different case and must
not silently replace this example's settings.

The existing `test/beddoesleishmanADG_tests.jl` uses separate high-precision
`data/ADG_intermediates.txt`: 48 reference times, 19 stations and 42 fields.
It includes separated flow, but its largest reference vortex contribution is
only about 0.106% of peak total normal force. The original regression defaults
to sound speed 343 m/s while its reference records 335 m/s.

Four vortex-related source divergences from OpenFAST v3.3.0 remain under review.
This documentation repair does not change those equations. Passing the
regression or running the example is not a general dynamic-stall validation.
See [Audit and Fix Plan](@ref) for quantitative evidence and questions for Adam.
