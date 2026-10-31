# Contributor changes and verification

Based on `fix/narrow-airfoilvector-getproperty` at
`15f75cd1800ed9d1b40b3312afb7a78a5c832faa`. Tested with Julia 1.13.1 on Windows.
The scientific audit and unresolved numerical questions are published in
[Audit and Fix Plan](docs/src/audit_and_fix_plan.md).

## Changes

- Prefer `solve_discrete`; `solve_indicial` remains the identical function object.
- Replace obsolete Getting Started calls with a complete BL-ADG example.
  Correct units, helper return values, paths, load ordering and coefficient rotation.
- Repair Øye initialization, separation return, attached-limit singularity,
  continuous geometry dispatch and coefficient extraction. Retain the existing
  Faber polynomial and lift-only convention. Unsupported normal-force mode
  raises an explicit error. Reactivate Øye tests and document limitations.
- Correct Gonzalez's load-count metadata to three; its aerodynamic recurrence
  and force equations are unchanged. Reference datasets and tolerances are unchanged.

## Executed checks

| Check | Result |
|---|---|
| Full active `test/runtests.jl` | 611 passed, 0 failed, 0 errored, 0 broken |
| New Øye numerical verification | 80 assertions included in the 611 |
| Literal Julia blocks from Getting Started, BL and Øye pages | 9 passed, 0 failed |
| Primary BL-ADG baseline | All 3,200 states and all saved coefficients reproduced exactly |
| Documenter 0.27.25 HTML build | Completed; 15 docstrings missing from manual reported as warnings |
| Earlier archived evidence | 295 SHA-256 hashes checked, none changed |
| Original master/fix checkouts | Clean; retained separately |

The pinned test runner retained all dependency versions/source trees from the
previous investigation, including SciMLBase 3.39.1, FLOWMath 0.3.3 and CurveFit
1.11.3. The new, separate example environment was independently instantiated;
it does not upgrade the research environments or tracked package manifests.
OpenFASTTools is pinned to `b16ff261855030b184868c1eafcc3ba8830b50b8`.

One initial documentation check failed: continuous versus discrete Øye
attachment differed by 0.178871 with unconstrained adaptive steps. `saveat`
does not constrain internal steps. Limiting `dtmax` to the input sample interval
resolved the check without changing model parameters or relaxing its 0.02 bound.
The initial failure log is retained in the separate local research workspace.
Visual review also caught an erroneous second degree-to-radian conversion in
the draft Øye example. The stored Faber polar already spans -pi to pi. Removing
that conversion restored plausible lift; a broad lift-range assertion now guards
against this mistake. Draft failure evidence was retained outside the package.

## Reproduce

Run from this checkout's root, with Julia 1.11+ (tested 1.13.1):

```sh
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples test/run_documentation_examples.jl
```

On this Windows installation add `--startup-file=no --pkgimages=no`. Set
`GKSwstype=100` for headless plots. The examples assert that the local source
is imported. Set `DOCUMENTATION_PLOT_DIR` to preserve generated plots after the
Julia process exits; otherwise the examples use temporary directories.

The full suite was executed with the preserved external `contributor-tests`
environment, invoking `test/runtests.jl` through a wrapper that recorded counts.
For a fresh test environment, activate a separate directory, develop this
checkout, add the pinned OpenFASTTools URL/revision above and the dependencies
in `test/Project.toml`, then include this checkout's `test/runtests.jl`.
This avoids changing the package's tracked project or manifest. Preserve the
resulting external Manifest to reproduce dependency resolution.

For an HTML build, use a separate environment containing this local checkout
and Documenter 0.27.25; set `DOCUMENTER_BUILD_ONLY=true` and include
`docs/make.jl`. This bypasses deployment. HTML output is ignored by Git.

## Scientific limits

Getting Started gives Cn RMSE **0.010582752184821022** against the supplied
100-row OpenFAST record. Cd/Cm discrepancies remain. Its low angles do not
exercise appreciable vortex loads. The unchanged separated-flow regression
still passes its existing tolerances, which do not establish vortex fidelity.
Øye analytic and interface tests establish numerical behavior, not experimental
dynamic-stall validation. The audit's four ADG source divergences require a
separate formulation decision and verification before numerical repairs.

No upstream merge, upstream push, pull request or external message is part of
this change. The audit is shared through the personal fork for review.
