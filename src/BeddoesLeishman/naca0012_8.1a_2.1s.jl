using DelimitedFiles, Polynomials, Plots, DifferentialEquations

include("../BeddoesLeishmanss.jl")

function U(t)
    M = 0.383
    a = 343.3
    return M*a
end

function Udot(t)
    return 0
end

function V(t)
    return 0
end

function Vdot(t)
    return 0
end

function alpha(t)
    c = 0.1
    M = 0.383
    a = 343.3
    shift = 2.1
    amp = 8.2
    k = 0.074

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function alphadot(t)
    c = 0.1
    M = 0.383
    a = 343.3
    shift = 2.1
    amp = 8.2
    k = 0.074

    v = M*a
    omega = k*2*v/c
    
    alfd = amp*omega*cos(omega*t)
    return alfd*(pi/180)
end

function alphaddot(t)
    c = 0.1
    M = 0.383
    a = 343.3
    shift = 2.1
    amp = 8.2
    k = 0.074

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
c = 0.1

### Polar Analysis
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/xf-n0012-il-1000000.csv", ','; skipstart=12)
polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989staticdata_fig1.csv", ',')


polar[:,1] = polar[:,1].*(pi/180) #Convert to radians

zeroalf, zeroidx = nearestto(polar[:,2], 0.0) #Nearest to zero lift
tenalf, tenidx = nearestto(polar[:,1], 10.0*pi/180) #Ten is typically before stall #*pi/180
linrng = zeroidx:tenidx

mxb = fit(polar[linrng,1], polar[linrng,2], 1)

dcldalpha = mxb.coeffs[2] #least squares fit 

alpha0 = roots(mxb)[1] #least squares fit   # NACA 0012 is a symmetric airfoil... so shouldn't the zero lift be at zero? 

# polarplt = plot(polar[:,1], polar[:,2], leg=false)
# plot!(mxb)
# display(polarplt)

maxcl, maxclidx = findmax(polar[:,2]) #critical lift coefficient (given by Larsen 2007), the beginning of LE separation. 
# println("Static Max Cl: ", maxcl)
alphav = polar[maxclidx, 1] 

maxcl = dcldalpha*(alphav-alpha0) #Projecting the inviscid lift to the static stall point
# println("Inviscid Max Cl: ", maxcl)

maxcl = 1.2 #Todo: I may need to rotate maxcl to Cn1

### Constants
M = 0.383
a = 343.3
v = M*a #Putting the speed here to calculate w. 


A = [0.3, 0.7, 1.0, 1.0]
b = [0.14, 0.53, 0.0, 0.0, 0.0] #I don't think these are currently used. 
S = [0.01, 0.03] #[2.75, 1.4].*(c/(2*v)) #TODO: I think this needs to be adjusted for each airfoil, so something matches the static data.
Tp = 1.7*(c/(1*v)) #1.7 
Tf = 3.0*(c/(1*v)) #3.0*(c/(1*v))
Tv = 6.0*(c/(1*v))
Tvl = 7.5*(c/(1*v))
a = 343.0 #TODO: I'm not entirely sure what this will do to the computations, but I think it'll have an effect. I'm not sure what to set it. 

### Initialize 
u0 = zeros(8)
u0[1] = 0.00
u0[2] = 0.00
u0[3] = 0.00
u0[4] = 0.00 #5e-6
u0[5] = 1.00 #Initial static value. 
u0[6] = 0.9 #Assuming the vortex starts attached
u0[7] = 0.00 
u0[8] = 0.00


p = [U, alpha, alphadot, c, dcldalpha, alphav, alpha0, maxcl, A, b, a, Tp, Tf, Tv]

dutest = zeros(8)
states_liftonly!(dutest, u0, p, 0.1)


tspan = (0.0, 1.0)

# dutest = zeros(8)
# states_liftonly!(dutest, u0, p, 0.1)

### Solve the system
prob = ODEProblem(states_liftonly!, u0, tspan, p)
sol = solve(prob) #;dtmax=0.0001

statesplt = plot(sol,linewidth=2,xaxis="t",title=["x1" "x2" "x3" "x4" "C\'n" "f\'\'" "tau" "Cvn"], leg=false,layout=(2,4))
# display(statespclt)

Cn, Cpn, Cdn, Ccn, Cina, Cinq, u, du, Cpn2, Ccndot, Cf, K = parsesolution_liftonly(sol, p)



alphavec = alpha.(sol.t)
alfavec = alphavec.*(180/pi)

exppolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/fig5rotated.csv", ',')

inviscid(alf) = dcldalpha.*(alf.-alpha0)
cli = inviscid(alphavec)

Cnplt = plot(alphavec, Cn, xaxis="α", yaxis="Cn", leg=false)
scatter!([alphavec[1]], [Cn[1]])

dispplot = plot(statesplt, Cnplt, layout=(2,1))
display(dispplot)

cycleplot = plot(alfavec[100:end], Cdn[100:end], lab="BL", leg=:topleft)
# plot!(alfavec, Cdn, lab="dynamic")
# plot!(alfavec, cli, lab="inviscid")
# plot!(aoavec.*(180/pi), clo, lab="inviscid")
scatter!(exppolar[:,1], exppolar[:,2], lab="experimental")
xlabel!("angle (degrees)")
ylabel!("Normal Coefficient")
xlims!((-10,15))
ylims!((-1.0,1.5))
display(cycleplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanss/leishman1989fig5/NACA0012_fig5_experimentalslope_082321a.png")

# Cdplot = plot(alpha.(sol.t).*(180/pi), Cdd)
# xlabel!("angle (degrees)")
# ylabel!("Coefficient of Pressure Drag")
# # display(Cdplot)
# # savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanlarsen/leishman1989fig5/NACA0012DRAG_fig5_experimentalslope_mintimestep0.0001.png")

# nothing