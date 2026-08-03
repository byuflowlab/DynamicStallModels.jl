# DynamicStallModels handoff — narrow the `getproperty(::AbstractVector{<:Airfoil})` override

**For:** an agent working in `~/.julia/dev/DynamicStallModels`
**Requested by:** Adam Cardoza · **Date:** 2026-07-31
**Context origin:** a WATT.jl performance change (store `Blade.airfoils` as a plain
`Vector{Airfoil}` instead of a `StructArray`) is blocked by a latent bug in this override.

## TL;DR

The override at `src/airfoils.jl:373-375`:

```julia
function Base.getproperty(obj::AbstractVector{<:Airfoil}, sym::Symbol)
    return getfield.(obj, sym)
end
```

is **too broad**. It claims `getproperty` for *every* `AbstractVector{<:Airfoil}` and *every*
symbol, and unconditionally broadcasts `getfield` over the elements. Under Julia ≥ 1.11,
`Array`'s own internals (`setindex!`, `collect`, `dataids`/`unalias`, broadcasting) access the
array's private fields **through `getproperty`** — e.g. `A.ref`, `A.size`. Those internal
accesses now hit this override, which tries to broadcast `getfield(airfoil, :ref)` over the
elements (a field that doesn't exist, on a possibly-undef vector) and throws.

**Net effect:** a plain `Vector{Airfoil}` is essentially unusable in Julia 1.11 — you cannot
even `setindex!` into one, nor `collect` one. This is why every caller currently must use a
`StructArray{Airfoil}` (StructArrays defines a *more specific* `getproperty(::StructArray, …)`
that wins dispatch, so the broad override never fires for a StructVector).

## Reproduction

```julia
using DynamicStallModels
const DS = DynamicStallModels
# build one real airfoil however you like, call it `af::DS.Airfoil`
v = Vector{DS.Airfoil}(undef, 3)
v[1] = af            # ERROR: UndefRefError
```

Minimal standalone reproduction (no DS needed — same override shape):

```julia
abstract type Foo end
struct Bar <: Foo; a::Int; b::Float64; end
Base.getproperty(obj::AbstractVector{<:Foo}, sym::Symbol) = getfield.(obj, sym)   # the bug

v = Vector{Foo}(undef, 3)
v[1] = Bar(1, 2.0)   # UndefRefError: setindex! -> v.ref -> getproperty -> broadcast over undef
```

### Observed stack traces (from the WATT test suite, Julia 1.11.9)

`setindex!` path:
```
UndefRefError: access to undefined reference
  [1] getindex                @ ./essentials.jl:917
  ...broadcast materialize...
  [8] getproperty(obj::Vector{Airfoil}, sym::Symbol)   @ DynamicStallModels ~/.julia/dev/DynamicStallModels/src/airfoils.jl:374
  [9] setindex!(A::Vector{Airfoil}, x::Airfoil{...}, i::Int64)   @ Base ./array.jl:987
```

`collect` path (building a `Vector{Airfoil}` from a `StructVector{Airfoil}`):
```
  [8] getproperty(obj::Vector{Airfoil}, sym::Symbol)   @ DynamicStallModels .../airfoils.jl:374
  [9] dataids(A::Vector{Airfoil})                      @ Base ./abstractarray.jl:1562
 [10] mightalias(A::Vector{Airfoil}, B::StructVector{Airfoil,…})   @ Base ./abstractarray.jl:1537
 [11] unalias ... copyto! ... collect
```

Both are Base reaching for `A.ref` / array identity via `getproperty` and getting hijacked.

## Root cause

`Base.getproperty(x, sym)` is the fallback for *all* dotted field access, and by default is
`getfield(x, sym)`. This override replaces that fallback for the whole `AbstractVector{<:Airfoil}`
type family, so it also intercepts the *structural* fields Julia's `Array` exposes to its own
runtime (`ref`, `size`). The intent of the override was only the **SoA convenience**: let user
code write `airfoilvector.polar` and get `[af.polar for af in airfoilvector]`. It should only do
that for symbols that are genuine `Airfoil` fields, and defer to normal `getfield` for anything
else.

`Airfoil` fields (from `fieldnames(Airfoil)`):
`(:model, :polar, :cl, :cd, :cm, :cn, :cc, :dcldalpha, :dcndalpha, :alpha0, :alphasep, :alphacut, :cutrad, :sfun)`

## Recommended fix (validated)

Narrow the override so it only field-broadcasts for real `Airfoil` fields and otherwise falls
back to the default `getfield` (i.e. normal Julia behavior for `Array.ref`, `Array.size`, etc.):

```julia
function Base.getproperty(obj::AbstractVector{<:Airfoil}, sym::Symbol)
    if sym in fieldnames(Airfoil)
        return getfield.(obj, sym)   # SoA convenience: broadcast a real Airfoil field
    else
        return getfield(obj, sym)    # everything else: Array internals (:ref, :size), default access
    end
end
```

- `fieldnames(Airfoil)` works on the parametric `Airfoil{…}` UnionAll and returns the 14 names
  above (verified). If you prefer to avoid the tuple `in` allocation on a hot path, hoist it to a
  `const` (`const _AIRFOIL_FIELDS = fieldnames(Airfoil)`) or use a `@generated`/`hasfield`-based
  check — but correctness-wise the simple version above is fine and `sym in fieldnames(Airfoil)`
  is constant-foldable.
- This is a **strict superset** of correct behavior: the `else` branch is exactly what an
  un-overridden `getproperty` would do, so nothing that currently works can regress. It only
  *adds* the ability for Base/`Array` internals (and hence plain `Vector{Airfoil}`) to function.

### Validated

The narrowed version was checked against the minimal reproduction: `setindex!`, `collect`, and
the SoA broadcast (`v.a`) all succeed:

```
NARROW override: setindex! OK
field-broadcast v.a = [1, 3, 5]
collect OK: 3
```

## Things to check / verify in the DS repo

1. **Grep for legitimate consumers of the SoA broadcast** (whole-vector `airfoils.<field>`):
   ```
   grep -rnE "\.(model|polar|cl|cd|cm|cn|cc|dcldalpha|dcndalpha|alpha0|alphasep|alphacut|cutrad|sfun)\b" src/ test/
   ```
   In `src/` there were **no** such internal uses (the override is used, if at all, by external
   callers / tests). Confirm the narrowed version still returns the broadcast for those symbols
   (it does).
2. **Run the DS test suite** — the fix should be behavior-neutral for all existing tests.
3. **Add a regression test** that a plain `Vector{Airfoil}` can be allocated-undef, `setindex!`-
   filled, iterated, and `collect`ed, and that `airfoilvector.polar` still returns the broadcast.
4. Consider whether the same issue affects any `AbstractMatrix{<:Airfoil}` / N-D usage; the fix
   generalizes if you also define/guard those, but the vector case is what WATT needs.

## Why WATT needs this

WATT is changing `Blade.airfoils` from a `StructArray{Airfoil}` (which reconstructs a fresh
immutable `Airfoil` on every `airfoils[i]` — ~34k reconstructions per AD march) to a plain
contiguous `Vector{Airfoil}` (O(1) pointer load). The WATT-side change (materialize once in the
`Blade` constructor; examples build a `Vector` directly) is written and ready, but it hits the
`UndefRefError` above the moment it tries to fill or store a plain `Vector{Airfoil}`. Once this
override is narrowed, the WATT change works and is bitwise-identical in results (the `Airfoil`
data is immutable and unchanged).

**When you're done:** report back that the override is narrowed and the DS tests pass; the WATT
session will then finish its Fix 1 (Vector storage) + Fix 2 (primal BEMT bracketing) + benchmark.
