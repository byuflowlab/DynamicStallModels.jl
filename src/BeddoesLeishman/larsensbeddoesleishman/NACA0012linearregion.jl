using DelimitedFiles
using Plots
using FLOWMath
using Polynomials

include("../BeddoesLeishmanLarsen.jl")

function b2w(b, v, c; a=343.3)
    #Not sure if this equation actually works. I got it from comparing wagner function approximations from Larsen and Leishman. 
    return b*(1-((v/a)^2))*2*v/c
end


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
tenalf, tenidx = nearestto(polar[:,1], 10.0*pi/180) #Ten is typically before stall
linrng = zeroidx:tenidx

mxb = fit(polar[linrng,1], polar[linrng,2], 1)

dcldalpha = mxb.coeffs[2]#least squares fit #Todo: There's an eta in Leishman 1989, that Larsen doesn't include. Where should that be accounted for. 

alpha0 = roots(mxb)[1] #least squares fit   # NACA 0012 is a symmetric airfoil... so shouldn't the zero lift be at zero? 

# polarplt = plot(polar[:,1], polar[:,2], leg=false)
# plot!(mxb)
# display(polarplt)

maxcl, maxclidx = findmax(polar[:,2]) #critical lift coefficient (given by Larsen 2007), the beginning of LE separation. 
# println("Static Max Cl: ", maxcl)
alphav = polar[maxclidx, 1] 

maxcl = dcldalpha*(alphav-alpha0) #Projecting the inviscid lift to the static stall point
# println("Inviscid Max Cl: ", maxcl)


### Constants
M = 0.383
a = 343.3
v = M*a #Putting the speed here to calculate w. 
w = ones(7) #Todo: I'm not really sure what to do about these constants. I suppose I could optimize them to match... but that seems like cheating. I could leave them the same and assume they will be close to that. 

w[1] = 0.125*2*v/c #I can rework this one if I have b1 #b2w(0.14, v, c) #Didn't do too much
w[2] = 0.375*2*v/c #I can rework this one if I have b2 #b2w(0.53, v, c)
w[3] = 0.275*2*v/c
w[4] = 0.075*2*v/c
w[5] = 2.5*2*v/c
w[6] = 2.5*2*v/c
w[7] = 0.4*2*v/c

A = [0.3, 0.7, 1.0, 1.0]
b = [0.14, 0.53] #I don't think these are currently used. 
S = [0.05, 0.05] #Todo: I think this needs to be adjusted for each airfoil, so something matches the static data. 
a = 343.0 #TODO: I'm not entirely sure what this will do to the computations, but I think it'll have an effect. I'm not sure what to set it. 

### Initialize 
u0 = zeros(8)
u0[1] = 0.0  
u0[2] = 0.0 
u0[3] = 0.00
u0[4] = 0.0 
u0[5] = 0.0
u0[6] = 0.75
u0[7] = 0.00 #Assuming the vortex starts attached
u0[8] = 0.00

states!, p, u0 = BeddoesLeishman(U, Udot, V, Vdot, alpha, alphadot, alphaddot, c, dcldalpha, alpha0, alphav, maxcl, A, S, w; u0=u0)
tspan = (0.0, 0.1)

### Solve the system
prob = ODEProblem(states!, u0, tspan, p)
sol = solve(prob;dtmax=0.0001)


# statesplt = plot(sol,linewidth=2,xaxis="t",title=["c1" "c2" "c3" "c4" "cl0d" "fd" "tau" "Clv"], leg=false,layout=(2,4))
# display(statesplt)

Cld, Cdd, Ccd = parsesolution(sol, p)
alphavec = alpha.(sol.t)

exppolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/fig5rotated.csv", ',')

inviscid(alf) = dcldalpha.*(alf.-alpha0)
cli = inviscid(alphavec)

cycleplot = plot(alpha.(sol.t).*(180/pi), Cld, lab="BLL", leg=:topleft)
# plot!(alpha.(sol.t).*(180/pi), cli, lab="inviscid")
# plot!(aoavec.*(180/pi), clo, lab="inviscid")
scatter!(exppolar[:,1], exppolar[:,2], lab="experimental")
xlabel!("angle (degrees)")
ylabel!("Coefficient of Lift")
# ylims!((0.0,1.5))
display(cycleplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanlarsen/leishman1989fig5/NACA0012_fig5_experimentalslope_mintimestep0.0001.png")

Cdplot = plot(alpha.(sol.t).*(180/pi), Cdd)
xlabel!("angle (degrees)")
ylabel!("Coefficient of Pressure Drag")
# display(Cdplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanlarsen/leishman1989fig5/NACA0012DRAG_fig5_experimentalslope_mintimestep0.0001.png")

nothing
