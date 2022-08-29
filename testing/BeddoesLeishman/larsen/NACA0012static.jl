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

function alphadot(t)
    return 0
end

function alphaddot(t)
    return 0
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

# ffunc = createf(polar)

zeroalf, zeroidx = nearestto(polar[:,2], 0.0) #Nearest to zero lift
tenalf, tenidx = nearestto(polar[:,1], 10.0*pi/180) #Ten is typically before stall
linrng = zeroidx:tenidx

mxb = fit(polar[linrng,1], polar[linrng,2], 1)

dcldalpha = mxb.coeffs[2]#least squares fit
#0.113*(180/pi) #leishman1989 state space paper
#mxb.coeffs[2]#least squares fit #Todo: There's an eta in Leishman 1989, that Larsen doesn't include. Where should that be accounted for. 

alpha0 = roots(mxb)[1] #least squares fit
# 0.17*(pi/180) #Leishman 1989 state space paper
# roots(mxb)[1] #least squares fit   # NACA 0012 is a symmetric airfoil... so shouldn't the zero lift be at zero? 

# polarplt = plot(polar[:,1], polar[:,2], leg=false)
# plot!(mxb)
# display(polarplt)

maxcl, maxclidx = findmax(polar[:,2]) #critical lift coefficient (given by Larsen 2007), the beginning of LE separation. 
# println("Static Max Cl: ", maxcl)
alphav = polar[maxclidx, 1] 
# 14.0*(pi/180) #Leishman 1989 state space paper
# polar[maxclidx, 1] 

maxcl = dcldalpha*(alphav-alpha0) #Projecting the inviscid lift to the static stall point


A = [0.3, 0.7, 1.0, 1.0]
b = [0.14, 0.53] #I don't think these are currently used. 
S = [0.025, 0.025] #Leishman 1989 A state space... 
# [2.75, 1.4] #Leishman 1989 state space paper
# [0.05, 0.05] #Leishman 1989 A state space... 
a = 343.3 

### Constants
M = 0.383
a = 343.3
v = M*a  

w = ones(7)  
w[1] = 0.125*2*v/c 
w[2] = 0.375*2*v/c 
w[3] = 0.275*2*v/c
w[4] = 0.075*2*v/c
w[5] = 2.5*2*v/c
w[6] = 2.5*2*v/c
w[7] = 0.4*2*v/c

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


aoavec = collect(-5:.1:20)
alphavec = aoavec.*(pi/180)
n = length(aoavec)
Clds = zeros(n)
Cdds = zeros(n)

for i = 1:n

    function alpha(t)
        return alphavec[i]
    end

    states!, p, uinit = BeddoesLeishman(U, Udot, V, Vdot, alpha,   alphadot, alphaddot, c, dcldalpha, alpha0, alphav, maxcl, A, S, w; u0=u0)
    tspan = (0.0, 0.3)

    ### Solve the system
    prob = ODEProblem(states!, uinit, tspan, p)
    sol = solve(prob;dtmax=0.0001)

    #Find the average state location
    Cld, Cdd, Ccd = parsesolution(sol, p)

    # if i==n #Plot the last solution. 
    #     alfavec = alpha.(sol.t).*(180/pi)
    #     plt = scatter(alfavec, Cld)
    #     hline!([Cld[end]])
    #     display(plt)
    # end

    Clds[i] = Cld[end]
    Cdds[i] = Cdd[end]

    #TODO: I could take the previous states and set them as the new initial state. 
end


inviscid(alf) = dcldalpha.*(alf.-alpha0)
cli = inviscid(alphavec)

Clplot = plot(aoavec, Clds, label="Dynamic Model", legend=:bottomright)
plot!(aoavec, cli, label="Insviscid")
scatter!(polar[:,1].*(180/pi), polar[:,2], label="Static")
xlabel!("Angle of Attack (Degrees)")
ylabel!("Coefficient of Lift")
display(Clplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanlarsen/leishman1989fig1/leishman1989ffunctionmatchedconstants.png")


nothing
