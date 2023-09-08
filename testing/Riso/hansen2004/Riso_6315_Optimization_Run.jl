using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations, SNOW

include("Riso_6315_Optimization.jl")

x_0 = [0.294, 0.331, 0.1, 0.3, 1.5, 6.0]

lx = [0.2, 0.05, 0.05, 0.05, 0.05, 0.05]

ux = [0.5, 0.5, 0.5, 0.5, 10.0, 10.0]

lg = [-Inf]

ug = [0.0]

ng = 1

ip_options = Dict("tol" => 1e-6)

options = Options(solver=IPOPT(ip_options), derivatives=CentralFD())

xopt, fopt, info = minimize(objective, x_0, ng, lx, ux, lg, ug, options)

### Optimized values
"""
A1 = 0.27469922114926987
A2 = 0.05433173899834567
b1 = 0.11224679387440323
b2 = 0.45235838971138226
Tp = 0.8173872939878122
Tf = 3.671051417177959

fopt = 0.02843954047511188
"""