# Beddoes-Leishman 
## Original and State Space Implementations
These implementations are still under construction. 


## AeroDyn Implementations
This model is also indicial form, so it will behave similarly to all the other indicial implementations. There are some key differences in this implementation from the original implementation, and they are noted in the theory tab. 

Since this same model is implemented in AeroDyn, we use OpenFASTsr.jl and some AeroDyn inputs and outputs to compare against. First we'll read in information about the airfoil. 

```julia
using DynamicStallModels, OpenFASTsr
of = OpenFASTsr
du21_a17 = of.read_airfoilinput("../../../data/DU21_A17.dat")
af = of.make_dsairfoil(du21_a17)
```

Now we create a vector of airfoils. Here we will only test a single airfoil, so we will store that lone airfoil in a vector. 
```julia
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af
```


Next we create an instance of the model that we want to simulate. We'll use Beddoes-Leishman, with Gonzalez's modifications (the third model). 
```julia
dsmodel = BeddoesLeishman(Indicial(), 1, airfoils, 3)
```

Now we need some inflow information. We'll use information that is passed to the unsteady module in OpenFAST. We need the inflow velocity and relative angle of attack for every time step of our simulation. 

```julia
using DelimitedFiles
fullout = readdlm("../../../data/aerodynout_fordynamicstall.csv", ',')
names = fullout[1,:]
names[1] = "Time"
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in 1:length(names)) 

tvec = outs["Time"]
nt = length(tvec)

Uvec = [sqrt(outs["AB1N011Vx"][i]^2 + outs["AB1N011Vy"][i]^2) for i in 1:nt] #m/s
aoavec = outs["AB1N011Alpha"].*(pi/180)
```



Before we can solve the time series, we need the chord length. Here we'll assume the length is 1m. Then we can solve. 

```julia
c = 1.0

states, loads = solve_indicial(dsmodel, [c], tvec, Uvec, aoavec)
```

Now we can visualize the outputs using `Plots`. 

```julia
Cn = loads[:,1] #The columns of loads are Cn, Cc, Cl, Cd, Cm

using Plots, LaTeXStrings

cnplt = plot(xaxis="Time (s)", yaxis=L"C_n")
plot!(tvec, Cn, lab="BL_AD")
plot!(tvec, outs["AB1N011Cn"], linestyle=:dash, lab="OpenFAST")
display(cnplt)
```
