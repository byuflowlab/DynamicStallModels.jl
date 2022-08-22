# DynamicStallModels.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](http://flow.byu.edu/DynamicStallModels.jl/)


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


The dynamic stall models are implemented in three ways: as a functional, an iterative, and in an indicial form. 

Functional implementations create a struct that returns the state rates. This allows for easy solution using DifferentialEquations.jl. These functional implementations are designed to be passed a single set of parameters, $p$. Out of this set of parameters, several are functions, such as the freestream velocity. 

Iterative implementations create a struct that returns the state rates as well, however they are not as easily solved by DifferentialEquations.jl. The input parameters are now all constant values, including the freestream velocity and angle of attack. This implementation is designed to be solved iteratively, meaning that for a given set of parameters, the states are updated for a single time step. The purpose of this type of implementation is for interfacing with other packages, specifically [Rotors.jl](https://github.com/byuflowlab/Rotors.jl). In the future, a function will be introduced that converts functions that describe environmental inputs into parameters, and then iterates through the solution of the model. 

The final implementation is an indicial formulation. Rather than providing a state space model to be solve, the model takes the environmental inputs and time step and calculates the states at the next time step. Several of the dynamic stall models included in this pakage were first developed in indicial form. 

Check out the documentation on how to get started and summaries of the implementation of the theories. 