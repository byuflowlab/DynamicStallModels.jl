using Plots, DifferentialEquations, FLOWMath, DelimitedFiles, Polynomials

function objective()
### Constants
M = 0.3
a = 343.3
v = M*a #Putting the speed here to calculate w. 

A = [0.3, 0.7, 1.0, 1.0]
b = [0.14, 0.53, 0.0, 0.0, 0.0]
S = [2.75, 1.4] 

T_p = 1.7
T_f = 3.0
T_v = 6.0
T_vl = 7.5

dCndalpha = 0.113*180/pi
alpha1 = 14.0*(pi/180)
alpha0 = 0.17*(pi/180)
Cn1 = 1.31

### Initialize 
u0 = zeros(8)
u0[1] = 0.00
u0[2] = 0.00
u0[3] = 0.00
u0[4] = 5e-6
u0[5] = 0.00
u0[6] = 1.00 #Assuming the vortex starts attached
u0[7] = 0.00 #Assuming the vortex starts attached
u0[8] = 0.00

tspan = (0.0, 0.5)

Ufun, Udotfun, alphafun, alphadotfun, alphaddotfun = prepenvironment(; c=chord, M=0.3, a=343.3, amp=10.0, shift=10.0, k=0.1)

conv1(t) = (c1/(2*Ufun(t)))

p = [Ufun, Udotfun, alphafun, alphadotfun, chord, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv1]

prob = ODEProblem(states!, u0, tspan, p)
sol = solve(prob;dtmax=0.0001)

t, u, du, Cnt, Cnp, Cdn, Cfc = parsesolution2(sol, p)

end