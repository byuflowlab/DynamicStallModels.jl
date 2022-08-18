# DynamicStallModels.jl

This package is a collection of different dynamic stall models. 

Dynamic stall models included in this package: 
- Beddoes-Leishman - State Space
- Beddoes-Leishman - Indicial
- Beddoes-Leishman - AeroDyn implementation
- Risø (Hansen 2004) - State Space
- Risø - Indicial
- Larsen
- Onera
- Oye


To get started, add the package. 
```julia
pkg> add https://github.com/byuflowlab/DynamicStallModels.jl.git
```

Then checkout the [Getting Started](@ref) page. 


!!! warning "Warning - Lack of Validation"
    An important thing to note is the lack of validation. The Beddoes-Leishman implementations are very far from validated and likely still have bugs. 
    The Risø models perform much better. For attached conditions, the models perform as expected compare excellently to published data (see the validation section of theory (that has yet to be created)). In stall conditions, when actual stall occurs, the models do not perform as expected. We are unsure if it has something to do with our implementation, or due to a lack of information provided in the published paper. We were unable to find the operating conditions that Hansen did his simulations at, and so we assumed them. We were able to get a decent match by optimizing the dynamic coefficients, but they were much different than the ones Hansen provided. 
    The Onera model has been validated (The validation needs to be added). 
