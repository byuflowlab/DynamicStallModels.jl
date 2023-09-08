using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, LaTeXStrings, DifferentialEquations, SNOW

include("Faber_0030_Optimization.jl")

x_0 = [0.197, 0.568, 0.05, 1.09, 0.97, 6.06]

lx = [0.02, 0.02, 0.02, 0.05, 0.05, 0.05]

ux = [0.5, 0.7, 0.7, 1.5, 10.0, 10.0]

lg = [-Inf]

ug = [0.0]

ng = 1

ip_options = Dict("tol" => 1e-8)

options = Options(solver=IPOPT(ip_options), derivatives=CentralFD())

xopt, fopt, info = minimize(objective, x_0, ng, lx, ux, lg, ug, options)

"""
Best Results:

A1 = 0.11901046494134122
A2 = 0.2964635220443628
b1 = 0.03833505254532422
b2 = 1.0129092704149298
Tp = 1.2267722555383287
Tf = 6.257533490189598

fopt = 0.04400000359837085
"""