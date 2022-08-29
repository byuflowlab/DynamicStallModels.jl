using Plots, DifferentialEquations, FLOWMath, DelimitedFiles, Polynomials

include("blss.jl")

function U(t)
    M = 0.3
    a = 343.3
    return M*a
end

function Udot(t)
    return 0
end

function alpha(t)
    c = .1
    M = 0.3
    a = 343.3
    shift = 10.0
    amp = 10.0
    k = 0.1

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function alphadot(t)
    c = .1
    M = 0.3
    a = 343.3
    shift = 10.0
    amp = 10.0
    k = 0.1

    v = M*a
    omega = k*2*v/c
    
    alfd = amp*omega*cos(omega*t)
    return alfd*(pi/180)
end

function alphaddot(t)
    c = .1
    M = 0.3
    a = 343.3
    shift = 10.0
    amp = 10.0
    k = 0.1

    v = M*a
    omega = k*2*v/c
    
    alfdd = -amp*(omega^2)*sin(omega*t)
    return alfdd*(pi/180)
end

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

### Geometry
c = .1
# k = 0.1
# M = 0.3
# a = 343.3
# f = 8.0
# c = k*M*a/(f*pi) #Todo: I'm not sure of this. This doens't seem like it should truly apply

### Polar Analysis
polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989staticdata_fig1.csv", ',')


polar[:,1] = polar[:,1].*(pi/180) #Convert to radians

zeroalf, zeroidx = nearestto(polar[:,2], 0.0) #Nearest to zero lift
tenalf, tenidx = nearestto(polar[:,1], 10.0*pi/180) #Ten is typically before stall #*pi/180
linrng = zeroidx:tenidx

mxb = fit(polar[linrng,1], polar[linrng,2], 1)

dCndalpha = mxb.coeffs[2] #least squares fit 

alpha0 = roots(mxb)[1] #least squares fit   # NACA 0012 is a symmetric airfoil... so shouldn't the zero lift be at zero? 

# polarplt = plot(polar[:,1], polar[:,2], leg=false)
# plot!(mxb)
# display(polarplt)

maxcn, maxcnidx = findmax(polar[:,2]) #critical lift coefficient (given by Larsen 2007), the beginning of LE separation. 
# println("Static Max Cl: ", maxcl)
alphav = polar[maxcnidx, 1] 

maxcn = dCndalpha*(alphav-alpha0) #Projecting the inviscid lift to the static stall point
Cn1, Cn1idx = nearestto(polar[:,2], maxcl)
alpha1 = polar[Cn1idx, 1]
# println("Inviscid Max Cl: ", maxcl)

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
a = 343.0

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

p = [U, Udot, alpha, alphadot, c, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a]

dutest = zeros(8)
states!(dutest, u0, p, 0.1)

### Solve the system
prob = ODEProblem(states!, u0, tspan, p)
sol = solve(prob;dtmax=0.0001) #;dtmax=0.0001

statesplt = plot(sol,linewidth=2,xaxis="t",title=["x1" "x2" "x3" "x4" "C\'n" "f\'\'" "tau" "Cvn"], leg=false,layout=(2,4))
display(statesplt)

t, u, du, Cnt, Cnp = parsesolution(sol, p)


alphavec = alpha.(t)
alfavec = alphavec.*(180/pi)

exppolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/leishman1989statespace/fig9/naca0012/experimental data.csv", ',')

modpolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/leishman1989statespace/fig9/naca0012/Model.csv", ',')

cycleplt = plot(legend=:topleft, xaxis="Angle of Attack (degrees)", yaxis="Normal Coefficient", title="NACA 0012")
scatter!(exppolar[:,1], exppolar[:,2], lab="Experimental")
scatter!(modpolar[:,1], modpolar[:,2], lab="Model - paper")
plot!(alfavec, Cnt, lab="Beddoes-Leishman")
display(cycleplt)



