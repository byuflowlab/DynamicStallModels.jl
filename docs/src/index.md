# DynamicStallModels.jl

This package is a collection of different dynamic stall models. 

The currently exercised interfaces are discrete AeroDyn/Gonzalez BL, continuous
and discrete lift-only Øye, continuous aggregate Risø with explicit constants,
and discrete NoModel static coefficients. Historical source files for additional
formulations remain in the repository; they are not all loaded or implemented
through the public solver interfaces. See [Audit repairs and compatibility](@ref).


To get started, add the package. 
```
pkg> add https://github.com/byuflowlab/DynamicStallModels.jl.git
```

Then checkout the [Getting Started](@ref) page. 


!!! warning "Warning - Lack of Validation"
    The examples and tests provide limited numerical verification and
    computational regression. They do not establish broad experimental
    dynamic-stall validation. The documented BL record stays at low incidence;
    historical separated-flow regression has weak vortex coverage and unresolved
    formulation/reference differences. Risø's attached-flow equation check does
    not validate separated stall or establish empirical time-constant calibration.
    No validation claim is made here for dormant Larsen/Onera/BL prototypes.


### References

Please check our References page to see the original authors of content presented herein. 


