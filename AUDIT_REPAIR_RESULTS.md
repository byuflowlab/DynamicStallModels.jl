# Audited defect repairs

Base: `acdf9bbc0425125c535bae76bfc24d2865276b95`. Repairs are on
`codex/audit-repairs`, separate from all research/audit baseline checkouts.
Tested on Windows with Julia 1.13.1. No upstream changes or publication were made.

## What changed

| Audit findings | Repair / scope |
|---|---|
| B01, B06 | Scalar and vector callbacks consistently use output-first ordering. Fixed high-level driver at the same time, state slicing and missing vector return. Added input-preservation and analytical-response tests. |
| B02, B03, B08 | Owned four-column radian polar, preserved supplied Cm, zero Cm only when omitted, consistent slope/stall metadata and updatable simple airfoils. Input matrices are preserved. |
| B04, B10 | Bound/validated recovery efficiency, handled the attached chordwise singularity, and replaced removed CurveFit.linear_fit with explicit ordinary least squares. No arbitrary slope fallback. |
| B05 | BL nonmutating update now uses its model argument and passes geometry correctly; compared with in-place and scalar paths. |
| B07 | Included the exact positive-stall angle in BLSP's exponential branch; both stall boundaries return .7 continuously. |
| B09 | Checked interpolation lengths, strict ordering, finite data and NaN queries. Typed interpolation fields. |
| B11 | Removed three undefined exports; exported actual make_airfoil/make_simpleairfoil helpers. No unsupported old constructor alias was invented. |
| B12 | Removed machine-specific subproject manifests, declared local/pinned sources, removed unused test ImplicitAD dependency, and repaired portable test/doc setup. Root package manifest retained. |
| B13 | OSP uses the scalar lift slope rather than calling it as a function. |
| B14 | Finite simultaneous zero-static/zero-inviscid lift limit in Risø Clfs. |
| B15 | Generic load wrapper returns the allocated three coefficients. |
| B16 | ADGSP lower-arity separation delegates to the slope-aware formula using the static normal slope. Dynamic ADG still supplies its own circulatory slope. |
| B17, B18 | Dormant legacy deficiencies/vortex helpers use decaying exponentials; reduced distance integrates 2U/c with trapezoidal quadrature, exact for linear speed histories. Separate from active ADG. |
| B19 | Standalone Larsen RHS explicitly receives A and theta with its other inputs; no undeclared globals. The incomplete old six-entry vector is rejected. This source remains experimental. |
| B20 | Finite geometry/environment, positive chord, positive subsonic ADG speed, positive ADG update time and correct state dimensions. Zero-speed model closure is explicitly unsupported. |
| B21 | Retired the unused private get_state_types helper rather than inventing an undocumented returnstatetype API. No repository callers existed. |
| B22, B23 | Risø residuals accept the airfoil/constants explicitly, use the correct rate signature and return dx-rates without mutating dx. Fixed the duplicate residual index. |

The aggregate continuous Risø wrapper additionally receives A/b/T explicitly,
uses each station's own states during extraction, and accepts numbers/functions
for environmental inputs. T is dimensional seconds; it is not silently converted
from unrelated nondimensional constants. The old introductory flat-plate inputs
are reproduced on the Risø example page, while Getting Started remains BL-ADG.

Unsupported BL variants now fail with an actionable error. Their incomplete
branches were removed from shared dispatch, not presented as working models.
NoModel remains a discrete static-polar model with complete common metadata.

## Executed verification

| Check | Observed result |
|---|---|
| New regressions before implementation edits | 2 passed, 28 failed, 7 errors; preserved failure log |
| Final native Pkg.test() on an isolated current-source snapshot | **726 passed, 0 failed, 0 errored** |
| Literal documentation code, freshly instantiated examples project | **12 passed** across Getting Started, BL-ADG, Øye and Risø |
| Risø attached-flow analytical comparison | Maximum Cl error **1.6268653091344731e-12** |
| ADG primary reproduction | Every saved state and coefficient identical to the preserved baseline |
| Original Gonzalez regression | All 456 existing assertions retained and passing; reference tolerances unchanged |
| HTML documentation build | Succeeded; 15 missing-manual-docstring warnings remain |
| Development environment | Fresh pinned resolution imports local DynamicStallModels and OpenFASTTools |

The final test source snapshot was hash-compared against the working repair
checkout. Passing tests cover specified properties and inputs; they do not prove
that every repository function or aerodynamic model is validated.

## Dependencies and reproducibility

Test/example/development projects use Julia 1.11+ sources to select the local
checkout and OpenFASTTools revision `b16ff261855030b184868c1eafcc3ba8830b50b8`.
GXBeam is pinned to tested version 0.1.2 in those environments: an unconstrained
development graph selected GXBeam 1.1.2 and ReverseDiff 1.18.7, which collided
during precompilation. No third-party source was patched. The earlier failed
resolution/import logs and generated unconstrained manifest are preserved.

The constructor repair was exercised with CurveFit 1.11.3 and 0.6.1 snapshots;
it no longer uses that package's removed linear_fit function. Historical research
environments were not upgraded. Generated subproject manifests are ignored;
retain local copies to reproduce a resolved dependency graph.

From this checkout's root:

```sh
julia --startup-file=no --pkgimages=no --project -e 'using Pkg; Pkg.test(;julia_args=`--startup-file=no --pkgimages=no`)'
julia --startup-file=no --pkgimages=no --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --pkgimages=no --project=examples test/run_documentation_examples.jl
```

Headless plots used GKSwstype=100. Documentation was built with
DOCUMENTER_BUILD_ONLY=true and `julia --project=docs docs/make.jl` after
instantiation. CI now checks supported test-environment Julia versions and runs
literal examples. The duplicate malformed documentation workflow was consolidated.
Cloud CI, other operating systems and Julia 1.11 were not executed locally.

## Remaining limits

ADG force/state recurrence equations and reference data are unchanged; the
previous vortex/OpenFAST equivalence questions remain open. The 726 assertions
are numerical/interface verification and existing computational regression,
not broad experimental dynamic-stall validation.

Historical development scripts with obsolete calls, missing optional packages
or hardcoded paths have not all been modernized. The full semantic/scientific
audit still has UNVERIFIED rows. These repairs address the confirmed findings,
without claiming that all legacy prototypes now work or that all possible bugs
have been removed. Earlier audit reports remain historical evidence.

Detailed local scripts, before/after error logs, native test logs, environment
snapshots and plots are retained in the separate `../repair_work` folder.
