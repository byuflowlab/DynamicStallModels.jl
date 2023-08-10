# Risø 

---

## Functional Solve Example

This version of the Risø method uses the functional solve. The functional solve uses an ODE solver like `DifferentialEquations.jl` to automaticall solve the state rate equations.

To begin, we prepare for a functional solve. In this example, we will be using the NACA 0030 polar from Faber's paper. This is a polar that was created using `webplotdigitizer.com` and was extended using extrapolation with `CCblade.jl`. `DelimitedFiles.jl` is then used to read in this file.

``` julia
using DynamicStallModels, DelimitedFiles
dsm = DynamicStallModels
file = "../../polars/Extended_NACA_0030.csv"
polar = readdlm(file , ',')
```

We then generate the constants that are used in this solve: such as chord length, mach number, and the inflow velocity.

``` julia
c = 0.55
M = 0.11
v = 343
Vrel = M*v
```

The model can now be generated. We set `Riso` as the dynamic stall model; `Functional()` as the solve method; and fill the vectors with the airfoil specific constants $A1$, $A2$, $b1$, $b2$, $T_p$, and $T_f$, respective from left to right.

Then the airfoil can be created by using the `dsm.make_airfoil` function and by passing is the polar, the model, and choosing the Risø separation function `dsm.RSP()`.

``` julia
dsmodel = Riso(Functional(), [0.11901046494134122, 0.2964635220443628] , [0.03833505254532422, 1.0129092704149298], [1.2267722555383287, 6.257533490189598]) 
af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.RSP())
```

An array can then be created (of type `Airfoil`) to hold the airfoil information, and then the single airfoil is then passed into the array.

``` julia
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af
```

The time range that the solve will be ran over can then be generated.

``` julia
tspan = (0.0, 2.0)
```

Then functions for the inflow velocity (`Uvector`), angle of attack (`alphavec`), and rate of change of angle of attack (`alphavecdot`) can be created.

``` julia
function Uvector(t)
    return Vrel
end

function alphavec(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.64
    k = 0.051

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function alphavecdot(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.64
    k = 0.051

    v = M*a
    omega = k*2*v/c

    alf = amp*omega*cos(omega*t)
    return alf*(pi/180)
end
```

The parameters can now be passed into a vector that will be used for the solve. The order for the vector must always be inflow velocity, rate of change of inflow velocity, angle of attack, and rate of change of the angle of attack (in that order). If multiple airfoils are be evaluated, the this pattern will just continue on: making the vector larger and larger. For this specific solve, the rate of change of the inflow velocity is zero, so `0.0` is passed into its position in the vector.

For the initial conditions, the order is $x_1$, $x_2$, $x_3$, and $x_4$. Passing in `0.0` for the initial conditions is perfectly fine, since the solution will converge to solutions anyways. Also, just like the parameter vector, if multiple airfoils are being evaluated, then the order just continues in this vector.
 
``` julia
parameters = [Uvector, 0.0, alphavec, alphavecdot]
x_initial = [0.0, 0.0, 0.0, 0.0]
```

The ODE problem can then be created using the `ODEProblem` function.

``` julia
prob = ODEProblem(airfoils, x_initial, tspan, parameters)
```

And the solution can be found using the `DifferentialEquations.solve` function.

``` julia
sol = DifferentialEquations.solve(prob, reltol=1e-8)
```

The `parsesolution` function can then be used to turn the matrix of state values from the solve function into aerodynamic coefficients and their corresponding angle of attacks.

``` julia
answer = parsesolution(dsmodel, airfoils, sol, parameters)
```

`Plots` can then be used to view the results of the solve.

``` julia
using Plots
plot(answer[1,:].*180/pi, answer[2,:], xlabel = "Angle of Attack (Degrees)", ylabel = "Cl", linewidth=3, label = "DSM (Risø)", color=:blue)
```

---

## Indicial Solve Example

The initial set up of this solve is the same at the `Functional()` solve (see above); however the first difference arises when creating the model. Instead of `Functional()`, use `Indicial()`.

``` julia
dsmodel = Riso(Indicial(), [0.11901046494134122, 0.2964635220443628] , [0.03833505254532422, 1.0129092704149298], [1.2267722555383287, 6.257533490189598]) 
af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.RSP())
```

An array of type `Airfoil` can then be created to store our airfoil in.

``` julia
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af
```

Another difference between this solve method and the previous one is the way that the time range is created. Instead of a simple tuple, a `range` function is used.

``` julia
tvec = range(0.0, 2.0, 1000)
```

`0.0` gives the starting time; `2.0` gives the ending time; and `1000` determines how many equally spaced values there are in between the starting and ending point (the starting and ending points account for two of these `1000` values).

The same process for generating the parameters in the `Functional()` solve can be used in the `Indicial()` solve.

``` julia
function Uvector(t)
    return Vrel
end

function alphavec(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.64
    k = 0.051

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function alphavecdot(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.64
    k = 0.051

    v = M*a
    omega = k*2*v/c

    alf = amp*omega*cos(omega*t)
    return alf*(pi/180)
end
```

Differently, though, these functions are passed into the solve not as their original functions but as vectors.

``` julia
avec = alphavec.(tvec)
uvec = Vrel.*ones(length(tvec))
avecdot = alphavecdot.(tvec)
```

The `solve_indicial` function can then be used to solve for the state values and loads.

``` julia
states, loads = solve_indicial(airfoils, tvec, uvec, avec; alphadotvec = avecdot)
```

`Plots` can then be used to view the results.

``` julia
using Plots
plot(avec.*180/pi, loads[:,1], linewidth=3, xlabel = "Angle of Attack (Degrees)", ylabel = "Cl", label = "DSM (Risø)", color=:blue)
```
