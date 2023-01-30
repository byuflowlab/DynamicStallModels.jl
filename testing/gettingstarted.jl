

k = 0.2 #reduced frequency
u = 1.0 #freestream velocity (m/s) 
c = 1.0 #Chord (m)
amplitude = 2
shift = 5

omega = k*2*u/c #Rate of oscillation

### Prepare inputs
U(t) = u #Freestream velocity as a function of time
Udot(t) = 0.0 #Derivative of freestream velocity as a function of time
alpha(t) = (amplitude*sin(omega*t) + shift)*(pi/180) #Angle of attack as a function of time (radians)
alphadot(t) = (amplitude*omega*cos(omega*t)*(pi/180)) #Pitching rate as a function of time (radians)

aoa = -pi:0.01:pi #Angle of Attack, (radians)
lift = 2*pi.*(aoa) #Coefficient of lift
drag = zero(aoa) #Coefficient of drag, left zero for this demonstration. 
polar = hcat(aoa, lift, drag)

A = [0.165, 0.335] #From the Hansen 2004 paper, for flat plate
b = [0.0455, 0.3000] # "" ""

Tp = 3.0 # "" "" 
Tf = 6.0 # "" "" 
T = [Tp, Tf]

using DynamicStallModels

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = airfoil(polar; A, b, T)

dsmodel = Riso(airfoils; detype=Functional())

using DifferentialEquations
x0 = zeros(4) #Hansen 2004's suggest initial conditions.
x0[3] = 1.0

p = [U, Udot, alpha, alphadot, c]

tspan = (0.0, 80.0)

prob = ODEProblem(dsmodel, x0, tspan, p)

sol = solve(prob, dtmax=0.1)

Cl, Cd, t =  parsesolution(dsmodel, sol, p)