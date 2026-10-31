# Audit and Fix Plan

## Scope and evidence

This contributor investigation reproduced the documentation, compared upstream
master `6cc8a23f0362a13e7479151068515354eac4298f` with Adam Cardoza's fix
`15f75cd1800ed9d1b40b3312afb7a78a5c832faa`, and examined Gonzalez separated-flow
references. Tests used Julia 1.13.1. OpenFASTTools was pinned to
`b16ff261855030b184868c1eafcc3ba8830b50b8`. No OpenFAST executable was run.

The master Getting Started page stopped at undefined `airfoil`; independent
probes also found obsolete tags and Risø wrappers accessing removed fields.
An external scalar adapter solved the original attached-flow Risø equations
and agreed with an analytic Cl solution to about 1.42e-12. The entire 803-row
solution CSV was identical on master and Adam's fix. His vector repair passed
16 targeted assertions and the 531-assertion active suite at the audited commit.
These counts describe the previous source snapshot, not all possible models.

The adapted DU21 Gonzalez example produced finite 100×32 states and 100×3
Cl/Cd/Cm loads. Correct Cn rotation gave full-record RMSE 0.0105828. Its largest
error occurred at 0.1 s; drag and moment did not agree tightly. Replacing chord
and sound speed with values from another repository case, plus a hypothetical
two-static-row startup, reduced Cn RMSE to 0.000370493, but did not establish
equivalent generating conditions or a repaired model. No reference rows were
discarded. The exact legacy CSV generator remains unknown.

## Vortex audit: unresolved numerical changes

The unchanged Gonzalez regression at the audited fix commit passed 456
assertions. At its DU40 node 5, reference peak Nv was 0.001604957 versus model
0.000163115; Nv RMSE 0.000508956 passed the test's absolute tolerance of 0.001.
The reference's strongest vortex contribution was only 0.106% of total peak Cn.
Close total loads therefore did not establish correct vortex dynamics.

Source comparisons use the Gonzalez branch of OpenFAST v3.3.0 at
[`5f3fb6ef74f48e75ca94000090737a41866fb264`](https://github.com/OpenFAST/openfast/blob/5f3fb6ef74f48e75ca94000090737a41866fb264/modules/aerodyn/src/UnsteadyAero.f90).
With F=(1+2sqrt(fpp))/3, four source-level divergences were identified:

| Mechanism | Audited Julia implementation | OpenFAST target |
|---|---|---|
| Vortex source | Ncirc*(1-F)^2 | Ncirc*(1-F²) |
| Clock/reset | Reset at 1+Stsh/Tvl under LESF | Gonzalez-specific gating involving Tvl+2(1-fpp)/Stsh |
| Decay modifier | Variable sigma3 | sigma3=1 for Gonzalez |
| Vortex moment | Positive xcp*(1-cos(pi*tauv/Tvl))*Nv | Negative corresponding term |

These are confirmed source divergences, not proven intentional variations.
The source expression also disagrees with the stored detailed C_V identity.
An **external frozen-state recurrence**, not a patched simulation, reduced Nv
RMSE to 0.000207417 when using the target source expression. It did not resolve
all discrepancies. The clock assertion in the existing test is disabled.
The original regression helper defaults to 343 m/s, while reference data record
335 m/s; it also recomputes Cn1 rather than using the file value directly.

This branch does **not** change those ADG equations, reference data or thresholds.
They need a decision on the intended formulation and separate verification.

## Changes in this contributor branch

1. Add `solve_discrete` as the preferred public name; retain `solve_indicial`
   as the identical function object, preserving existing calls.
2. Repair the shared continuous ODE dispatch/geometry interface required by
   Øye, and accept concrete Airfoil vectors in the discrete solver.
3. Repair Øye's attached-flow singularity and numeric separation return,
   align continuous/discrete time constants and extraction, and reactivate
   numerical verification tests. Keep the original lift-only Cd=Cm=0 convention;
   reject the previously empty normal-force mode instead of inventing closure.
4. Convert Getting Started to BL-ADG, correct the ADG and Øye examples/theory,
   and test the literal documentation code in a separate example environment.
5. Correct Gonzalez's load-count metadata to three. No ADG force/state formula
   is changed. See `CONTRIBUTOR_VALIDATION.md` in the repository root for the
   new branch's actual check results and reproduction commands.

## Questions for Adam and next work

- Which original OpenFAST revision, UA variant, polar files and startup schedule
  produced `aerodynout_fordynamicstall.csv`? Should Gonzalez use another dataset?
- Are the four vortex departures intentional, and which equations/sign
  conventions should govern a proposed numerical repair?
- Which sound speed/Cn1 convention should define regression equivalence?
- Is there a provenance-complete pitching reference with appreciable vortex
  loads, initial states, full-precision inputs and intended error tolerances?
- Should the next validation priority be attachment/separation recovery or
  vortex onset/convection/shedding?

The present evidence supports executable cases and limited computational
regression. It does not establish broad experimental dynamic-stall validation.
