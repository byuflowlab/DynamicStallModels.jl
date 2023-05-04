# Øye 

---

## Indicial Solve Example
---
This version of the Øye method utilizes the indicial solve. An indicial solve analytically solves the state rate equation and solves for new state values with a time step. 

To begin, we prepare for our indicial solve. In this example, a NACA 0015 polar from Faber's paper is read in using the package `DelimitedFiles`. Additionally, the model that we will be using for the solve is set up in this block. The Øye method is chosen for the dynamic stall model; an indicial solve is picked; the number `1` indicates that the coefficient of lift is being evaluated; the number `2` shows that the Faber/Larsen use of Hermite interpolation for the fully separated lift is used; and the `4.0` value is the constant used in the time scale parameter.

```julia
using DynamicStallModels, DelimitedFiles
dsm = DynamicStallModels
polar_0015 = readdlm("NACA_0015_Faber.csv" , ',')
dsmodel = Oye(Indicial(), 1, 2, 4.0)
```
Next, the chord length, mach number, and velocity are generated.
```julia
c = 0.55
M = 0.11
a = 343.0
Vrel = M*a 
```
The airfoil function is now created. In order to do this, the polar, the model, and the chord length are passed into the `make_airfoil` function. The `sfun` value at the end indicates that we are using the Larsen separation point function for this solve. 

An array is then generated to hold the information of all the airfoils that we want to evaluate. In this case, only one airfoil is being tested, so this array is of length 1 and filled with our single airfoil function.

``` julia
af = dsm.make_airfoil(polar_0015, dsmodel, c; sfun=dsm.LSP())
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af
```
Here, the range of time values is set up, and the angle of attack function with respect to time is created.
```julia
tvec = range(0, 1.0, 1000)
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
With the previous values prepared, a vector of the angle of attack at each time value is created. Additionally, the inflow velocity vector is made. Since the inflow velocity for this problem is constant, every value in this vector will be the same.
```julia
alphavec = alpha.(tvec)
Uvec = Vrel.*ones(length(tvec))
```
The state values and loads can be found using the `solve_indicial` function.
``` julia
states, loads = solve_indicial(airfoils, tvec, Uvec, alphavec)
```
This block of code changes the load values to either be with respect to lift or normal force coefficients.
``` julia
cn = loads[:,1]
if dsmodel.cflag == 2
    cn_static = af.cn.(alphavec)
else
    cn_static = af.cl.(alphavec)
end
```
Now, the dynamic loads or state values can be viewed using `Plots`
``` julia
using Plots, LaTeXStrings
stateplt = plot(tvec, states[:,1], leg=false, xaxis="time (s)", yaxis="f")
cyclecnplt = plot(xaxis="Angle of Attack (deg)", yaxis=L"C_L", leg=:bottomright)
plot!(alphavec.*(180/pi), cn, lab="DSM")
plot!(alphavec.*(180/pi), cn_static, lab="Static")
display(cyclecnplt) 
```

