using DynamicStallModels, OpenFASTTools, FLOWMath, Plots, Plots.PlotMeasures, DelimitedFiles, LaTeXStrings, Revise, Statistics
using BenchmarkTools, Traceur
#=
Time the BLADG functions. 

#Todo: 
- Rewrite the code to be designed for timing. 
- Time the initial speed. 
- Check for type stability
- Check memory management


Adam Cardoza 1/15/24
=#

#Todo: alpha (generic function with 5 methods). I'm not sure that there should be 5 methods. 

of = OpenFASTTools
DSM = DynamicStallModels

path = dirname(@__FILE__)
cd(path)


 

### Create airfoil
c = 3.542
idx = 1

du21_a17 = of.read_airfoilinput("../../../data/airfoils/DU40_A17.dat")


# af = Airfoil(polar, clfit, cdfit, cmfit, dcndalpha, alpha0, alphasep, A, b, T, S, xcp)
af = of.make_dsairfoil(du21_a17, c)
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af






### Create simulation data
tvec = 0.:0.001:10.
nt = length(tvec)
Uvec = 10.0.*ones(nt)
aoavec = (5.0 .+ 3.0.*sin.(2*pi.*tvec/2.0)).*(pi/180)



### Solve
states, loads = solve_indicial(airfoils, tvec, Uvec, aoavec)

### Before code optimization: 
@btime solve_indicial(airfoils, tvec, Uvec, aoavec) 
### 193.910 ms (5,436,334 allocations: 104.71 MiB) 
# @code_warntype(solve_indicial(airfoils, tvec, Uvec, aoavec)) #Outputs seem unstable. 
# @trace solve_indicial(airfoils, tvec, Uvec, aoavec) -> 
### After fixing the type concreteness of Airfoil, we got 15.413 ms (170,092 allocations: 10.46 MiB)... so an order of magnitude speed up and memory decrease.... which is bueno. 

### Looking at the ADG function. 
airfoil = af
oldstates = states[end-1, :]
newstates = zeros(32)
y = [Uvec[end], 0.0, aoavec[end], 0.0]
deltat = tvec[end]-tvec[end-1]
DSM.update_states_ADG!(airfoil, oldstates, newstates, y, deltat)

# @btime DSM.update_states_ADG!(airfoil, oldstates, newstates, y, deltat) ### 441.333 ns (0 allocations: 0 bytes)
# @code_warntype(DSM.update_states_ADG!(airfoil, oldstates, newstates, y, deltat))




### Testing the rotating the loads
# Cn = zero(tvec)
# Cc = zero(tvec)

# for i in eachindex(tvec)
#     Cn[i], Cc[i] = DSM.rotate_load(loads[i,1], loads[i,2], states[i, 1])
# end 





nothing