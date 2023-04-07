# Separation Function Notes

**Important Info:**
The separation function setup section starts on line 251 of [airfoils.jl](../../src/airfoils.jl) with the line `alpha0, _ = brent(cl, -0.25, 0.25)` (*Note: I don't know what this does*)

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

The if statement checks if sfun is of the ADSP type, and then runs the line given above. This is only necessary if ADSP is getting called as it is the only separation point function that needs some prework (I believe)

### To Do's

- [x] linear_fit doesn't actually return anything, line 278
- [ ] Write a test file like beddoesleishmanADG_tests.jl for the seperationpoint functions
  - [ ] take a look at test_seperationpoint.jl 
  - [x] Maybe make a variable that contains all the flags needed to clarify which method throughout the whole process so I only need to change 1 thing **Ans:** Adam doesn't want this and would like to keep the flexibility of choosing whatever you want and mixing/matching
  - [ ] Look at Weston's Documentation to better understand

## Notes for the different separation point functions

*These start around line 392, but before getting there it goes through the ADSP function*  
**ADSP:** If ADSP (at line 74) (AeroDyn Separation Point) is called it directly calls the `reverse_separationpointcalculation_ADO()` function (at line 518)  
**Larsen:** starts around line 639

- [ ] make sure the separation point functions in airfoils are correct
- [x] update oye update_states_oye_hansen the `fst = ...` line (about 189) to call the separation point function that I need `separationpoint(airfoil,....)` it doesn't need to specify type as it has multiple dispatch

## Notes for Oye.jl

*Questions:*  

- What is detype::DEType ?  
  - Ans: the type of model it is, Functional(), Iterative(), or Indicial(). We are just using Indicial()
- cflag::TI - what does this, and similar lines do? around line 25, does it say the data type for the struct?
  - Ans: Yesish, it names that field and gives it a datatype. But in this case, rather than a specific datatype, it just has to be the same datatype as the other fields. 
- How does it know which separation point function to use? or is that decided earlier?

*Notes:*  

- the `cflag` field says what says whether to apply the separation delay to lift (1) or normal force (2)  
- the `version` field is for Hansen(1), Faber (2) of the Oye model

## Notes for the whole path from OyeComparer.jl

1. Setup: the packages are imported and named, the polar is read in, and the airfoil is created.
   
   1. `dsm.airfoil(Polar; A = 8.0, sfun=dsm.LSP())` this is the line that generates the airfoil. `sfun=dsm.LSP()` tells it which separation point function to use
   
   2. | SP type | Input                                                                                                               | Explanation                                                                                                                                                                                                                                                                                                                                                                                                                           |
      | ------- | ------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
      | SP      | `dsm.SP=(1,1)` -> the inputs will become an akima spline ffit                                                       | Separation point function created by Adam Cardoza. It calls the `reverse_separationpointcalculation` function and then puts an Akima spline on it                                                                                                                                                                                                                                                                                     |
      | BLSP    | `dsm.BLSP(1)`  -> becomes a vector of floats holding S constants that are a best fit for the separation point curve | Beddoes-Leishman original separation point function. It uses a series of if statements to check what region we are in, and then assigns a value  - [  ] figure out if it is looking for and finding the region it thinks it is supposed ot                                                                                                                                                                                            |
      | ADSP    |                                                                                                                     |                                                                                                                                                                                                                                                                                                                                                                                                                                       |
      | ADGSP   |                                                                                                                     |                                                                                                                                                                                                                                                                                                                                                                                                                                       |
      | RSP     |                                                                                                                     | Same as HSP                                                                                                                                                                                                                                                                                                                                                                                                                           |
      | OSP     |                                                                                                                     | Parabolic fit                                                                                                                                                                                                                                                                                                                                                                                                                         |
      | HSP     |                                                                                                                     | Same as RSP                                                                                                                                                                                                                                                                                                                                                                                                                           |
      | LSP     | `dsm.LSP()` -> no inputs                                                                                            | Larsen separation point function from 2007 paper, and repeated in Oye from Faber 2018. Line 639 in airfoils.jl It calls `cl_fullysep_faber` to calculate the cl (or cn) values and then uses the same calculation as Oye, Larsen, and Faber $fst = (C_l^{st} - C_l^{fs}) / (C_l^{inv} - C_l^{fs})$ . st = static, fs = fully separated, inv = inviscid. It then clips values larger than 1 or smaller than 0... Hermite Interpolation |
      
      2. It than makes the Oye model struct `dsmodel = Oye(Indicial(), 1, airfoils,2,2)` this calls the solver type, how many airfoils, the airfoils, a flag for if the delay is applied to $C_l$ (1) or $C_n$ (2), and then a flag for solving with the Hansen (1) or Faber (2) model
      
      3. `solve_indicial(dsmodel....)` takes us to the solver
      
      4. in the indicial solver it calls the `initialize` function which goes to Oye and uses the flags and if statements to run the initial separation point calculation
      
      5. In solve it then continues and gets to `update_environment!(dsmodel,....)` - [ ] I don't really know why or what this does
      
      6. It eventually makes a round about way to Oye.jl lines 208 and uses the faber/larsen method to calculate the separation point. *Note:* that tau is hardcoded as `tau = A*c/U` 
      
      7. It then makes its way to `getcoefficient_indicial_faber` function in Oye.jl line 295, here It calculates Cn using Hansen's equation for Cn_inv and final (returned) Cn, and Faber for the fully separated Cn, *can they play all together like that?* 
      
      8. It then returns the states and loads and goes all the way back to my OyeComparer.jl file

# Will Do's

- [ ] Make the Larsen separation point functions and  related Oye solver fully functional and multiple dispatch
  
  - [ ] Fix separation point functions
    
    - [ ] Make sure I can call the Larsen function and get into it where I need to
      
      - [ ] this is linked with multiple dispatch down below, I can implement it I believe
      - [ ] Delete HSP
      - [ ] Add a flag to all separation point functions for Cl or Cn, default for ADSP and ADGSP and Riso/Hansen (RSP/HSP) -> needs to be Cn. Larsen's paper uses Cl
      - [ ] Create the change A function
    
    - [x] make sure I am using the larsen separation point and not hansen's when it expects it
      
      - [x] find the difference between the two, one overshoots, and one under?
      - [x] Ans: Larsen and Faber both want the fully separated alphasep, so in deep stall, in the Larsen vertol case, that is 32 deg
      - [x] Use update airfoil function to make sure alpha0 is the same as in Larsen's paper
      - [x] try 1/tau and just straight tau for the omega3 = 0.07-> Solved?: my conversion was wrong
      - [ ] Add and check my static Cl plot and webplot digitize his polar/static plot (they are the same)
      - [x] We want to be able to mismatch functions, un-hardcode Hansen's function so it is a function call
    
    - [x] I implemented `dsmodel.version == 3` as being Larsen and Adam already has it as BeddoesLeishman, fix all of them to `dsmodel.version == 4`
    
    - [x] Make sure all the equations in the Larsen function are correct
    
    - [x] Check the update states (step 6) equations and Tau
  
  - [ ] Create a separation point test file (do this first and match figure 4-c)
    
    - [x] Write out the questions I want to answer, what the answers should be, and where/when to find them
    
    - [x] create the file and testset
    
    - [ ] make it automated
    
    - [x] test against figure 4-c (Larsen Paper) separation point curve
  
  - [x] Fix Oye solver functions (check with Weston on this)
    
    - [x] check the equations in step 7 above
  
  - [x] Make multiple dispatch
    
    - [x] Add the `separationpoint` function call in the 1 or 2 places necessary, comment out the if statements
    
    - [x] check how the Oye solvers are called and fix/update them (check with Weston during this also)

### Grad School Project Ideas

- Dr. Ning's Lab
  
  - Modeling propeller and turbine wakes with the vortex particle method (Eduardo may have already done this?)
  
  - Talk with Ryan and Taylor, Tyler?
    
    - Ryan mentioned the BEM rotor interaction problem
  
  - Judd -> might have some ideas, but different than what I have done before -> panel method = panel to velocity to pressures to lift