# Seperation Function Notes 
**Important Info:**
The seperation function setup section starts on line 251 of [airfoils.jl](../../src/airfoils.jl) with the line `alpha0, _ = brent(cl, -0.25, 0.25)` (*Note: I don't know what this does*)

**Setup Line by line analysis:**  
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
Puts the minimum alpha value and the maximum alpha value in a vector. The alpha values appear to still be in radians  

```Julia
middlepolar = polar[minclidx:maxclidx,:]
```  
This creates a smaller polar just between our min and max cl values, and contains both the alpha and cl values. *Note: Issues were happening here as for the vertol polar the minclidx is after the maxclidx which causes problems, so maybe it is the findmin and findmax functions/lines that we need to fix.*  

```Julia
_, cl0idx = nearestto(middlepolar[:,2], 0.0)
```  
This calls the nearestto function and tries to find the Cl value in the smaller middle polar (between min and max Cl) closest to 0.0. The nearestto function subtracts the value we are looking for from the polar, takes the absolute value and then finds the minimum (both value and index) in the new vector. It takes the index value and uses that to find the minimum Cl in the original middle polar, the function returns that value and index (for the middlepolar).  

```Julia
alpha50 = middlepolar[end,1]*0.25
```  
This finds the alpha value at the end of the middle polar (so Clmax) and multiplies it by .25  

```Julia
_, alf50idx = nearestto(middlepolar[:,1], alpha50)
```  
This uses the nearestto function discussed above to find the index of the alpha value that is closest to alpha50, or the alpha value that is 1/4 of the Cl max alpha  

```Julia
    _, dcldalpha = linear_fit(middlepolar[cl0idx:alf50idx,1], middlepolar[cl0idx:alf50idx,2])
```  
This pulls the alpha and Cl values in from the minimum Cl alpha (in middle polar) up to alpha50 (1/4 of the Cl max alpha). The linear fit function then fits a line through the set of points, and returns the slope and y-intercept, but just the slope is saved as dclalpha  

```Julia
if isnan(dcldalpha)
    dcldalpha=2*pi #flat plate slope
    @warn("dcldalpha returned NaN")
end
```  
This sets the dcldalpha (or slope) to 2*pi, the standard flat plate slope if something goes wrong with the linear_fit, *which is happening*  

```Julia
sfun = ADSP(alphavec, cnvec, ccvec, alpha0, alphasep, dcldalpha, eta)
```  
The if statement checks if sfun is of the ADSP type, and then runs the line given above. This...  

- [x] linear_fit doesn't actually return anything, line 278
- [ ] Write a test file like beddoesleishmanADG_tests.jl for the seperationpoint functions
  - [ ] take a look at test_seperationpoint.jl 
  - [ ] Maybe make a variable that contains all the flags needed to clarify which method throughout the whole process so I only need to change 1 thing
  - [ ] Look at Weston's Documentation to better understand


## Notes for the different seperation point functions  
*These start around line 392, but before getting there it goes through the ADSP function*  
**ADSP:** If ADSP is called it directly calls the `reverse_separationpointcalculation_ADO()` function  
**Larsen:** starts around line 639


## Notes for Oye.jl  
*Questions:*  
- What is detype::DEType ?  
  - Ans: the type of model it is, Functional(), Iterative(), or Indicial(). We are just using Indicial()
- cflag::TI - what does this, and similar lines do? around line 25, does it say the data type for the struct?
  - Ans: Yesish, it names that field and gives it a datatype. But in this case, rather than a specific datatype, it just has to be the same datatype as the other fields. 
- How does it know which seperation point function to use? or is that decided earlier?

*Notes:*  
- the `cflag` field says what says whether to apply the seperation delay to lift (1) or normal force (2)  
- the `version` field is for Hansen(1), Faber (2) of the Oye model