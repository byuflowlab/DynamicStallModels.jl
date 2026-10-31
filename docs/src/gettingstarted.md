# Getting Started

This example runs the **Beddoes–Leishman AeroDyn/Gonzalez (BL-ADG)** model for
one DU21 airfoil using sampled inflow from the repository. It uses the discrete
recurrence, not an ODE solver. For continuous integration see the [Øye](@ref)
example. The previous Risø introduction used interfaces that no longer matched
the current Airfoil representation.

## Environment and data

Clone this repository and run commands **from its root directory**, which
contains `Project.toml` and `data/`. Use the separate example environment:

```sh
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples
```

The example Project uses Julia 1.11+ `[sources]` entries to select this local
checkout and pin OpenFASTTools. It declares plotting/data dependencies rather
than relying on your global environment. Instantiation creates a local example
Manifest; keep it to repeat your run. Do not activate `docs/` or `test/` as a
substitute. The development checks described below were run with Julia 1.13.1.
On Windows systems that block native package images, use `--pkgimages=no`.

Copy the following Julia blocks in order into that session. The data are read
only; plots are saved in a new temporary directory.

```julia
using DynamicStallModels, OpenFASTTools, DelimitedFiles
repo = pwd()
@assert isfile(joinpath(repo, "data", "DU21_A17.dat")) "Run from the repository root."
@assert realpath(pathof(DynamicStallModels)) == realpath(joinpath(repo, "src", "DynamicStallModels.jl"))

input = OpenFASTTools.read_airfoilinput(joinpath(repo, "data", "DU21_A17.dat"))
af, xcp = OpenFASTTools.make_dsairfoil(input)
airfoils = [af]
@assert af.model.detype isa Discrete
@assert af.model.version == 3
```

`make_dsairfoil` returns **two values**, the Airfoil and nondimensional vortex
centre parameter `xcp`. It also creates `af.model`; a separate
`BeddoesLeishman(Indicial(), ...)` call is neither necessary nor supported.
Julia model version 3 denotes Gonzalez; OpenFAST's corresponding switch is
UA_Mod=2. With this pinned helper the defaults are Akima polar interpolation,
sound speed 343 m/s and filter coefficient 0.5. `xcp` comes from the airfoil file.

## Prescribed inflow

The CSV has a header, units row, then 100 numeric rows. Its first header contains
a UTF-8 byte-order mark; normalizing it to `Time` avoids a mismatched key.
Node `AB1N011` supplies velocities in m/s and angle in degrees.

```julia
raw = readdlm(joinpath(repo, "data", "aerodynout_fordynamicstall.csv"), ',')
names = String.(raw[1, :])
names[1] = "Time"
data = Float64.(raw[3:end, :])
outs = Dict(names[i] => data[:, i] for i in eachindex(names))
t = outs["Time"]
U = sqrt.(outs["AB1N011Vx"].^2 .+ outs["AB1N011Vy"].^2)
alpha = deg2rad.(outs["AB1N011Alpha"])
chord = 1.0 # metres: the legacy documented example's choice
```

`t` spans 0–9.9 s at 0.1 s intervals. Angles passed to the model are **radians**.
The solver initializes from the first input; subsequent steps hold the next
sample's input constant over the step. The example's 1 m chord is not established
as the chord used to generate the reference CSV.

## Solve and extract loads

```julia
states, loads = solve_discrete(airfoils, t, U, alpha;
                               cvec=[chord], xcpvec=[xcp])
@assert size(states) == (length(t), 32)
@assert size(loads) == (length(t), 3)
@assert all(isfinite, states) && all(isfinite, loads)
Cl, Cd, Cm = loads[:, 1], loads[:, 2], loads[:, 3]
Cn = Cl .* cos.(alpha) .+ Cd .* sin.(alpha)
Ct = Cl .* sin.(alpha) .- Cd .* cos.(alpha)
```

The **three load columns are Cl, Cd, Cm**, per airfoil. `Cn` and `Ct` are total
normal/tangential-to-chord coefficients reconstructed by rotation; Cl is not Cn.
Cm uses the supplied quarter-chord, positive-nose-up convention. The 32 states
include lags and flags; they are not 32 continuously integrated ODE states.
`solve_indicial` remains an exact compatibility alias for `solve_discrete`.

## Compare and plot

```julia
reference_Cn = outs["AB1N011Cn"]
cn_rmse = sqrt(sum(abs2, Cn .- reference_Cn) / length(t))
@show cn_rmse

using Plots
fig = plot(layout=(3, 1), size=(900, 900))
for (i, label, calculated) in ((1, "Cn", Cn), (2, "Cd", Cd), (3, "Cm", Cm))
    plot!(fig[i], t, calculated; label="BL-ADG", ylabel=label, xlabel="Time (s)")
    plot!(fig[i], t, outs["AB1N011" * label]; label="Supplied reference", linestyle=:dash)
end
plot_path = joinpath(mktempdir(), "bl_adg_comparison.png")
savefig(fig, plot_path)
@assert isfile(plot_path)
@show plot_path
```

The audited primary run gave Cn RMSE approximately **0.01058**, with the largest
error near 0.1 s. It did **not** reproduce all coefficients tightly: Cd and Cm
disagree visibly. Do not remove startup rows, substitute Cl for Cn, or tune
parameters just to make the plot agree.

This record stays near 3–4° and does not activate meaningful dynamic stall.
Successful finite output verifies execution, not vortex physics. The exact
legacy CSV generator/UA variant/initial states remain unknown. See
[Beddoes–Leishman](@ref) and the [Audit and Fix Plan](@ref) for the separated-flow
regression, source divergences and unresolved scientific questions.
