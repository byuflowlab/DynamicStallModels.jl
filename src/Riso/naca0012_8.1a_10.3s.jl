using DelimitedFiles
using Plots
using Polynomials
using FLOWMath

include("../Riso.jl")


function U(t)
    M = 0.379
    a = 343.0
    return M*a
end

function Udot(t)
    return 0
end

# function V(t)
#     return 0
# end

# function Vdot(t)
#     return 0
# end

function alpha(t)
    c = 0.1
    M = 0.379
    a = 343.0
    shift = 10.3
    amp = 8.1
    k = 0.075

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function alphadot(t)
    c = 0.1
    M = 0.379
    a = 343.0
    shift = 10.3
    amp = 8.1
    k = 0.075

    v = M*a
    omega = k*2*v/c
    
    alfd = amp*omega*cos(omega*t)
    return alfd*(pi/180)
end

function alphaddot(t)
    c = 0.1
    M = 0.379
    a = 343.0
    shift = 10.3
    amp = 8.1
    k = 0.075

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




function alpha34(t) #Assuming that the aoa at the three-quarters point is the same as the aoa for now, but it really shouldn't be. 
    return alpha(t)
end

### Geometry
c = 0.1

### Polar Analysis
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/xf-n0012-il-1000000.csv", ','; skipstart=12)
polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989staticdata_fig1.csv", ',')


polar[:,1] = polar[:,1].*(pi/180) #Convert to radians

liftfit = Akima(polar[:,1], polar[:,2])

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
maxcl = 1.2
# println("Inviscid Max Cl: ", maxcl)


### Constants
M = 0.379
a = 343.0
V = M*a #Putting the speed here to calculate w. 

#Constants
A = [0.165, 0.335]
b = [0.0455, 0.3]
Tp = 1/0.4125
Tf = 1/0.0875


#Initialize 
x0 = zeros(4)
p = [c, A, b, Tp, Tf, dcldalpha, liftfit, U, Udot, alpha, alphadot, alpha0]
tspan = (0.0, 0.1)

prob = ODEProblem(states!,x0,tspan,p)
sol = solve(prob)

statesplt = plot(sol,linewidth=2,xaxis="t",label=["x1" "x2" "x3" "x4"],layout=(4,1))
display(statesplt)


Cld, f = parsesolution(sol, p)
alphavec = alpha.(sol.t)
Cls = liftfit.(alphavec)
alfavec = alphavec.*(180/pi)

clplt = plot(sol.t, Cld, lab="dynamic", leg=:bottomright)
plot!(sol.t, Cls, lab="static") #TODO: note that the static data only goes to around 15 degrees and so as the aoa oscillates, it is oscillating out of the region of provided data. 
xlabel!("time (s)")
ylabel!("Coefficient of Lift")
# display(clplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/riso/Cld_stepaoa.png")

fplt = plot(sol.t, f, leg=false)
xlabel!("time (s)")
ylabel!("Separation factor")
display(fplt)

cycleplt = plot(legend=false)
plot!(alfavec, Cld)
xlabel!("Angle of attack")
ylabel!("Coefficient of Lift")
display(cycleplt)

nothing