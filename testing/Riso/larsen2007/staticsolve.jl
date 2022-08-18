using DifferentialEquations, Plots, DelimitedFiles, FLOWMath, Roots, NLsolve

include("../riso.jl")


#=
Try solving the steady state solution of the Riso method using NLsolve. 
=#


# function static_solve(fun, x0, p, t, lowbounds, upbounds)
#     dx0 = zero(x0)

#     solveme! = function(resids, x)
#         fun(resids, dx0, x, p, t) #Todo: I'm running into typing issues.
#     end


#     mcpsolve(solveme!, lowbounds, upbounds, x0, autodiff = :forward)
# end





U(t) = 60
Udot(t) = 0.0

alpha(t) = 20*pi/180

alphadot(t) = 0.0



#Enviroment
V = 60

#Geometry
c = 1.5
polar = readdlm("../../../polars/extendedVertol 23010-1.58.dat", ',')
liftfit = Akima(polar[:,1], polar[:,2])
dragfit = Akima(polar[:,1], polar[:,3])
dcldalpha = 2*pi*1.05
alpha0 = 0.014 # -0.019
Cd0 = dragfit(alpha0)

alphas = [find_seperation_alpha(liftfit, dcldalpha, alpha0)...]
afm = alphas[2]
afp = alphas[1]

#Constants
A = [0.165, 0.335] 
b = [0.0455, 0.3]
Tp = 1/0.4125
Tf = 1/0.0875


#Initialize 
x0 = zeros(4)
x0[4] = 0.98
p = [U, Udot, alpha, alphadot, c, A[1], A[2], b[1], b[2], dcldalpha, alpha0, Tp, Tf, liftfit, afm, afp]

prob = SteadyStateProblem{true}(riso_ode!, x0, p)
sol = solve(prob)

# lb = [-Inf, -Inf, -Inf, 0.0]
# ub = [Inf, Inf, Inf, 1.0]

# result = static_solve(static_riso_dae!, x0, p, 300.0, lb, ub)


Cl, Cd = riso_coefficients(sol.u, U(0), alpha(0), alphadot(0), c, liftfit, dragfit, dcldalpha, dragfit(alpha0), alpha0, A[1], A[2], afm, afp)




nothing