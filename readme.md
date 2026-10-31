# DynamicStallModels.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](http://flow.byu.edu/DynamicStallModels.jl/)
[![Build Status](https://github.com/byuflowlab/DynamicStallModels.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/byuflowlab/DynamicStallModels.jl/actions/workflows/CI.yml?query=branch%3Amaster)


This package is a collection of different dynamic stall models.

Current execution paths in this checkout:

- Beddoes–Leishman AeroDyn/Gonzalez: discrete recurrence, 32 stored entries,
  three load columns (Cl, Cd, Cm).
- Øye: continuous or discrete attachment lag, lift only.
- Risø: continuous aggregate wrapper with explicit A, b and dimensional T.
- NoModel: discrete static-polar evaluation.

Other BL variants and continuous BL are unfinished and rejected by the shared
solver interfaces. Historical Larsen/Onera, BL state-space and other prototype
files exist, but are not loaded as supported package entry points. Their presence
does not imply validation. See [repair notes](docs/src/repair_notes.md) for limits.


Use [Getting Started](docs/src/gettingstarted.md) with the isolated example
environment. It selects the local checkout and pins the external helper.
Julia 1.11+ is required by the separate test/example project source settings.
`solve_discrete` is the sampled-input driver; `solve_indicial` remains an alias.
The scalar/vector callback contract is output first: `af(new, old, y, p, dt)`.

Passing tests verifies the exercised computations, not broad experimental
dynamic-stall fidelity. In particular, the low-angle example does not validate
vortex dynamics, and the separated-flow reference comparison has unresolved
formulation/configuration questions. See the [audit](docs/src/audit_and_fix_plan.md).

The repository also contains older functional, iterative and indicial prototypes.

Continuous interfaces return state rates for an ODE solver. Environmental
entries can be functions of time or numbers; the documented wrapper determines
how geometry and constants are supplied.

Discrete interfaces update states for sampled environmental inputs and a time
increment. `solve_discrete` provides this driver for the supported per-Airfoil
models, including the indicial recurrence used by BL-ADG.

Older scripts can use obsolete Functional/Indicial tags, missing helper context
or nonportable paths. They are historical development material; use the current
documented examples and [repair results](AUDIT_REPAIR_RESULTS.md) for tested paths.

Check out the documentation on how to get started and summaries of the implementation of the theories.
