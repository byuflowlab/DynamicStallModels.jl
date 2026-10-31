# Audit repairs and compatibility

These repairs address reproducible software and numerical defects. They do not
establish broad experimental validation or resolve the separate ADG vortex audit.

## State and coefficient interfaces

`af(new, old, environment, geometry, dt)` and
`airfoils(new, old, offsets, environment, geometry, dt)` are **output first**.
The discrete solver now uses that same ordering. The input state is preserved;
nonmutating calls allocate and return the result. Generic load extraction returns
all three coefficients instead of a mutation routine's incidental scalar return.

ADG's supported path is Discrete() version 3. Other BL variants remain unfinished
and now raise errors. Its mathematical domain requires positive chord, finite
inputs, positive time increments and 0<U<a. Zero-speed closure is not invented;
choose an appropriate physical treatment outside the model. These checks are not
evidence of accuracy throughout the entire subsonic range.

## Polars and compatibility

The convenience constructors copy input data, normalize angles to radians and
store alpha/Cl/Cd/Cm columns. Three-column inputs gain zero moment; four-column
inputs retain their moment. `radians=false` applies to the polar's angle column;
explicit cut angles remain radians. The lift slope is fitted with ordinary
unweighted least squares, independent of CurveFit's removed linear_fit function.
No 2pi slope is silently substituted for a failed fit.

The default ADSP recovery efficiency is the model's eta when present, otherwise
the same unit-efficiency default used by BeddoesLeishman; an explicit eta keyword
overrides it. The zero-incidence chordwise fit uses the existing attached limit.
No reference datasets or aerodynamic parameters were tuned.

Undefined `airfoil`, `simpleairfoil` and `riso` exports were removed. Use
`make_airfoil`, `make_simpleairfoil` and `Riso`. The unused private get_state_types
helper was removed rather than assigning it an undocumented numerical meaning.
`solve_indicial` still aliases `solve_discrete` exactly.

## Risø and experimental legacy code

The aggregate continuous Riso wrapper requires explicit A, b and T; T is in
seconds. It supports numeric/function-valued environments and separate station
states. Residuals receive the airfoil/constants explicitly and preserve dx.
See [Risø](@ref) for the former introductory case through the repaired public API.

The dormant indicial BL helpers now use decaying exponentials for homogeneous
deficiencies/vortex response, and trapezoidal integration of 2U/c for reduced
distance. These changes are separate from active ADG, whose recurrence equations
are untouched. The legacy source remains outside the package's include graph.

The standalone Larsen RHS takes a NamedTuple containing `CL0`, `CL0dot`, `alpha`,
`alphadot`, `alphav`, `w`, `A` and `theta`. `theta(t)` must supply the intended
static separation angle; amplitudes/closure are no longer undeclared globals.
The old incomplete six-entry vector raises an actionable error. This source
remains experimental and is not newly advertised as a validated public model.

## Environments and checks

Test, development and documentation projects select the local checkout through
Julia 1.11+ sources. OpenFASTTools is pinned to an immutable revision. Their
machine-specific tracked manifests were removed; locally generated manifests
are ignored and should be retained by the investigator for repeat runs.
The root package manifest and older research evidence are preserved.
GXBeam, a dependency of the pinned OpenFASTTools helper, is constrained to the
tested 0.1.2 in test/example/development environments. A fresh unconstrained
development resolution selected GXBeam 1.1.2 with ReverseDiff 1.18.7, which
failed precompilation by overwriting the same method. No third-party source was
patched and historical environment versions were not upgraded.

From the checkout root:

```sh
julia --project -e 'using Pkg; Pkg.test()'
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples test/run_documentation_examples.jl
```

On the tested Windows installation native package images were disabled with
`--pkgimages=no`; headless plotting used `GKSwstype=100`.
CI targets Julia 1.11/1.13 on Linux/Windows x64. The duplicate documentation
workflow was consolidated, fixing its path-filter/undefined-module defects.
Official documentation deployment remains restricted to upstream master pushes;
local verification uses DOCUMENTER_BUILD_ONLY=true.

See AUDIT_REPAIR_RESULTS.md in the root for executed results and remaining scope.
Historical development scripts with obsolete APIs, missing optional dependencies
or hardcoded paths are not all modernized by these repairs. Original audit
evidence remains unchanged, rather than being rewritten to show passing results.
