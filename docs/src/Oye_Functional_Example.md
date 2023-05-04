# Øye 
---
## Functional Solve Example

This version of the Øye method uses the functional solve. The functional solve uses an ODE solver like `DifferentialEquations.jl` to automaticall solve the state rate equations.

To begin, we prepare for the functional solve. In this example, a NACA 0015 polar from Faber's paper is read in using the package DelimitedFiles. Additionally, the model that we will be using for the solve is set up in this block. The Øye method is chosen for the dynamic stall model; the functional solve is picked; the number `1` indicates that the coefficient of lift is being evaluated; the number `2` shows that the Faber/Larsen use of Hermite interpolation for the fully separated lift is used; and the `4.0` value is the constant used in the time scale parameter.
```julia
using DynamicStallModels, DelimitedFiles
dsm = DynamicStallModels
polar_0015 = readdlm("NACA_0015_Faber.csv" , ',')
dsmodel = Oye(Functional(), 1, 2, 4.0)
```
The chord length, mach number, and inflow velocity are then prepared.
```julia
c = 0.55
M = 0.11
a = 343.0
Vrel = M*a 
```
The airfoil function is now created. In order to do this, the polar, the model, and the chord length are passed into the `make_airfoil` function. The `sfun` value at the end indicates that we are using the Larsen separation point function for this solve.

An array is then generated to hold the information of all the airfoils that we want to evaluate. In this case, only one airfoil is being tested, so this array is of length 1 and filled with our single airfoil function.
```julia
af = dsm.make_airfoil(polar_0015, dsmodel, c; sfun=dsm.LSP())
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af
```
A touple containing the span of time that we want to evaluate over is now generated.
```julia
tspan = (0, 2.0)
```
Two functions with respect to time are created, as well. `Uvector` gives the inflow velocity and `alpha` gives the angle of attack of the airfoil.
```julia
function Uvector(t)
    return Vrel
end
function alpha(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.0
    k = 0.051

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end
```
The parameters and initial condition for this solve can now be generated. The order for the parameter vector is imperative. The order for the inputted values should be the inflow velocity, the rate of change of the inflow velocity, the angle of attack, and the rate of change of the angle of attack. If multiple airfoils are being evaluated, then this pattern will just continue onwards through the vector's indices. Since this example only solves a single airfoil, only four inputs are required into the parameters vector. Additionally, the rate of change of the angle of attack and inflow velocity are not required for the Øye method, so they are just inputted in as `0.0`.

Additionally, any initial condition between `0.0` and `1.0` is good for the Øye method. The solution will converge if a value is picked in this range. Since we are only evaluating a single airfoil, there is only one initial condition required for this method.
``` julia
parameters = [Uvector, 0.0, alpha, 0.0]
x_initial = [0.8]
```
The ODE problem is now created using the `ODEProblem` function.
```julia
prob = ODEProblem(airfoils, x_initial, tspan, parameters)
```
The problem from the previous block of code can now be passed into the `solve` function of `DifferentialEquations.jl`. The keyword `reltol` adds more control to the step size in the solve and allows for more precise solutions.
```julia
sol = DifferentialEquations.solve(prob, reltol=1e-8)
```
With the solution found, the `parsesolution` function can be used to find the dynamic lift or normal force coefficients. This function outputs a matrix that contains alternating rows of the angle of attacks and their corresponding coefficient values for each airfoil.
```julia
answer = parsesolution(dsmodel, af, sol, parameters)
```
Using `Plots` allows us to view the results of our functional solve.
```julia
using Plots, LaTeXStrings
cyclecnplt = plot(xaxis="Angle of Attack (deg)", yaxis=L"C_L", leg=:bottomright)
plot(answer[1,:].*180/pi, answer[2,:])
```