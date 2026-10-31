# API Reference

Airfoils hold their model and static coefficient fits. Chord and xcp are
supplied to the solver, rather than stored as dynamic Airfoil fields.

## Airfoils

```@docs
DynamicStallModels.Airfoil
DynamicStallModels.make_airfoil
DynamicStallModels.make_simpleairfoil
DynamicStallModels.update_airfoil
```

## Solving

```@docs
solve_discrete
```

`solve_indicial` is the same function object as `solve_discrete` and remains
available for compatibility. New examples use the latter name.
For continuous Øye models use
`ODEProblem(airfoils,x0,tspan,y; cvec, xcpvec)` and an ODE solver, as shown in
[Øye](@ref). The environment has four entries per airfoil; geometry is separate.

## Beddoes–Leishman API

```@docs
BeddoesLeishman
```

The supported example uses the AeroDyn/Gonzalez discrete version. Read the
[Beddoes–Leishman](@ref) limitations before interpreting results as validation.

## Øye API

```@docs
Oye
```

## Risø API

```@docs
Riso
```

The continuous Risø aggregate wrapper requires explicit `A`, `b` and dimensional
`T` constants. See [Risø](@ref) for the attached-flow reproduction. It is separate
from the per-Airfoil discrete solver. Larsen/Onera helpers
that are not loaded by the module are not advertised as runnable entry points.
