# Seperation Function Notes 
**Important Info:**
The seperation function section starts on line 251 of [airfoils.jl](../../src/airfoils.jl) with the line `alpha0, _ = brent(cl, -0.25, 0.25)` (*Note: I don't know what this does*)

**Line by line analysis:**  
```Julia
alpha0, _ = brent(cl, -0.25, 0.25)
```  
*I am not sure yet* something about finding a root?  

```julia
_, minclidx = findmin(polar[:,2])  
_, maxclidx = findmax(polar[:,2])
```
Searches for a minimum (or maximum) cl value from the polar and returns the index  

```Julia
alphasep = [polar[minclidx, 1], polar[maxclidx,1]]
```  
Puts the alpha values between the min and max cl in a vector together