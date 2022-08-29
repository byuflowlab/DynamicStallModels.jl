using DelimitedFiles, Polynomials, Plots, DifferentialEquations

include("../BeddoesLeishmanss.jl")

function U(t)
    M = 0.379
    a = 343.0
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
maxcl = 1.2
# println("Inviscid Max Cl: ", maxcl)


### Constants
M = 0.379
a = 343.0
v = M*a #Putting the speed here to calculate w. 

A = [0.3, 0.7, 1.0, 1.0]
b = [0.14, 0.53, 0.0, 0.0, 0.0] #I don't think these are currently used. 
S = [0.01, 0.03] #[2.75, 1.4].*(pi/180) #Todo: I think this needs to be adjusted for each airfoil, so something matches the static data. 
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

#TODO: I may need to rotate maxcl to Cn1

p = [U, alpha, alphadot, c, dcldalpha, alphav, alpha0, maxcl, A, b, a, Tp, Tf, Tv, Tvl]


tspan = (0.0, 0.05)

# dutest = zeros(8)
# states_liftonly!(dutest, u0, p, 0.1)

### Solve the system
prob = ODEProblem(states_liftonly!, u0, tspan, p)
sol = solve(prob;dtmax=0.0001) #;dtmax=0.0001

statesplt = plot(sol,linewidth=2,xaxis="t",title=["x1" "x2" "x3" "x4" "C\'n" "f\'\'" "tau" "Cvn"], leg=false,layout=(2,4))
display(statesplt)

Cdn, Cpn, Cdd, Cdd2, Cdd3, Cdd4, Cdd5, Cdd6, Cdd7, Cdd8, Cdd9, Cdd10, Cdd11, Cdd12, Cdd13, Cdd14, Cdd15, Cdd16, Cdd17, Cdd18, Cdd19, u, du, Cc, Cfc, Cfc2, Cfc3, Cfc4, Cfc5, Cfc6, Cfc7, Cfc8, Cfc9, Cfc10, Cfc11, Cfc12, Cfc13, Cfc14, Cfc15, Cfc16, Cfc17, Cfc18, Cfc19 = parsesolution(sol, p)

u1plt = plot(sol.t, u[:,1], leg=false, title="x1")
u2plt = plot(sol.t, u[:,2], leg=false, title="x2")
u3plt = plot(sol.t, u[:,3], leg=false, title="x3")
u4plt = plot(sol.t, u[:,4], leg=false, title="x4")
u5plt = plot(sol.t, u[:,5], leg=false, title="C\'n")
u6plt = plot(sol.t, u[:,6], leg=false, title="f\'\'")
u7plt = plot(sol.t, u[:,7], leg=false, title="τ")
u8plt = plot(sol.t, u[:,8], leg=false, title="Cvn")
dynstatesplt = plot(u5plt, u6plt, u7plt, u8plt, layout=(4,1))
# display(dynstatesplt)


alphavec = alpha.(sol.t)
alfavec = alphavec.*(180/pi)

exppolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989fig8Cn.csv", ',')

mat = hcat(alfavec, Cdn)
# writedlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/outputs/leishman1989/naca0012_8.1_10.3s_Cdn_BLss.csv", mat, ',')

inviscid(alf) = dcldalpha.*(alf.-alpha0)
cli = inviscid(alphavec)

Cnpplt = plot(sol.t, u[:,5], xaxis="time (s)", yaxis="Cn\'", leg=false)
dCnpplt = plot(sol.t, du[:,5], xaxis="time (s)", yaxis="dCn\'", leg=false)
CnpdCnpplt = plot(Cnpplt, dCnpplt, layout=(2,1))
# display(CnpdCnpplt)

cycleplot = plot(alfavec[100:end], Cdn[100:end], lab="BL", leg=:topleft)
# plot!(alfavec, cli, lab="inviscid")
# plot!(aoavec.*(180/pi), clo, lab="inviscid")
# plot!(alfavec, Cdn, lab="Cdn")
plot!(alfavec, Cpn, lab="Cpn")
scatter!(exppolar[:,1], exppolar[:,2], lab="experimental")
xlabel!("angle (degrees)")
ylabel!("Normal Coefficient")
# xlims!((-10,15))
# ylims!((-1.0,1.5))
display(cycleplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanss/leishman1989fig8/NACA0012_fig8_experimentalslope_082321a.png")

cdexpdata = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989fig8Cd.csv", ',')

Cdplot = plot(leg=:topleft)
plot!(alfavec, Cdd, lab="Cdd") #Good match on the lower half of the returning path, but not on the outgoing path, and doesn't do well at all on the upper half. 
scatter!(cdexpdata[:,1], cdexpdata[:,2], lab="experimental") 
# plot!(alfavec, Cdd2, lab="Cdd2") #Same as below, but slightly higher
# plot!(alfavec, Cdd3, lab="Cdd3") #Correct max magnitude, wrong slope. 
# plot!(alfavec, Cdd4, lab="Cdd4") #Goes all over the place. 
# plot!(alfavec, Cdd5, lab="Cdd5") #Pretty negative. 
# plot!(alfavec, Cdd6, lab="Cn") #The lift times sin(α) #Cdd6 
# plot!(alfavec, Cdd7, lab="Cdd7") #Correct slope for the first half, but doesn't pick up at higher angles of attack. 
# plot!(alfavec, Cdd8, lab="Cdd8") #Larger and below experimental
# plot!(alfavec, Cdd9, lab="Cdd9") #Lower than Cdd
# plot!(alfavec, Cdd11, lab="Cdd11") #Good match on the returning curve, but way below for the outgoing curve.
# plot!(alfavec, Cdd12, lab="Cdd12") #Pretty much the same as Cdd.
# plot!(alfavec, Cdd10, lab="Cdd10") #The lowest when it comes to Cc*sin(α)
# plot!(alfavec, Cdd13, lab="Cdd13")
# # plot!(alfavec, Cdd14, lab="Cdd14") #Matches well on the lower half of the curve, both in and out, but not at the higher angles of attack.... Which makes me think that it's missing something with f. 
# # plot!(alfavec, Cdd15, lab="Cdd15") #lift*sin(alphaE)
# plot!(alfavec, Cdd16, lab="Cdd16")
# plot!(alfavec, Cdd17, lab="Cdd17")
# plot!(alfavec, Cdd18, lab="Cdd18")
plot!(alfavec, Cdd19, lab="Cdd19")
plot!(alfavec, Cfc19, lab="Cfc19")
xlims!((0,20))
ylims!((-0.1, 0.6))
xlabel!("angle (degrees)")
ylabel!("Coefficient of Pressure Drag")
# display(Cdplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanss/leishman1989fig8/NACA0012DRAG_fig8_experimentalslope_082421b.png")



### Trying to fit the Cdc data so it matches the experimental Cdd
aa = 1/300
bb = 0.55
cc = 0.01
eq(x) = aa*(x^2) + bb*x + cc
eqout = zeros(length(sol.t))

for i=1:length(sol.t)
    alf = alpha(sol.t[i])
    eqq = eq(alf)
    eqout[i] = Cdn[i]*sin(alf)-eqq*cos(alf)
end


expplt = plot(leg=:topleft)
scatter!(cdexpdata[:,1], cdexpdata[:,2], lab="experimental")
plot!(alfavec, eqout, lab="$aa x^2 + $bb x + $cc")
# xlims!((0,20))
# ylims!((-0.1, 0.6))
xlabel!("angle (degrees)")
ylabel!("Coefficient of Pressure Drag")
# display(expplt)


eqvec = eq.(alphavec)
Ccplt = scatter(alfavec, eqvec, lab="backworked CC", leg=:topleft)
# plot!(alfavec, Cc, lab="Cc") #above the line with that same weird peak at upstroke stall. 
# plot!(alfavec, Cfc, lab="Cfc") #Above the mark
# plot!(alfavec, Cfc2, lab="Cfc2") #to much angle, and the "cross" is at a later angle of attack. 
# plot!(alfavec, Cfc3, lab="Cfc3") #Straddles the line, but has a large increase on the upstroke. 
# plot!(alfavec, Cfc4, lab="Cfc4") #Good width, wrong angle
# plot!(alfavec, Cfc5, lab="Cfc5") #Good width, wrong angle, similar to 4, but lower. 
# plot!(alfavec, Cfc6, lab="Cfc6") #Much to high
# plot!(alfavec, Cfc7, lab="Cfc7") #Below at low angles, and above at high angles. It has an interesting peak around stall. 
# plot!(alfavec, Cfc8, lab="Cfc8") #Good width, wrong angle, similar to, but lower than 5
# plot!(alfavec, Cfc9, lab="Cfc9") #Too tall, but goes above and below.
# plot!(alfavec, Cfc10, lab="Cfc10") #Below for most then peaks above in upstroke. 
# plot!(alfavec, Cfc11, lab="Cfc11") #Too tall, but goes above and below. Similar to 9
# plot!(alfavec, Cfc12, lab="Cfc12") #Good width, wrong angle, similar to 8
# plot!(alfavec, Cfc13, lab="Cfc13") #Good width, wrong angle, basically the same as 12
# plot!(alfavec, Cfc14, lab="Cfc14") #Much too high, same as Cfc6
plot!(alfavec, Cfc15, lab="Cfc15") #Too tall but goes above and below
plot!(alfavec, Cfc16, lab="Cfc16") #The same 15
plot!(alfavec, Cfc17, lab="Cfc17")
plot!(alfavec, Cfc18, lab="Cfc18")
plot!(alfavec, Cfc19, lab="Cfc19")
xlabel!("Angle of Attack (degrees)")
ylabel!("Chordwise Coefficient of Force")
# display(Ccplt)
nothing