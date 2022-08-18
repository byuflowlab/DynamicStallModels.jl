function Heavi(t)
    if t>0
        return 1
    else
        return 0
    end
end

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

nothing

using Plots, DifferentialEquations, FLOWMath, DelimitedFiles, Polynomials
nothing

function prepenvironment(; c=0.1, M=0.3, a=343.3, amp=10.0, shift=10.0, k=0.1)
    v = M*a
    omega = k*2*v/c
    
    function U(t)
        return M*a
    end

    function Udot(t)
        return 0
    end

    function alpha(t)
        alf = shift + amp*sin(omega*t)
        return alf*(pi/180)
    end

    function alphadot(t)  
        alfd = amp*omega*cos(omega*t)
        return alfd*(pi/180)
    end

    function alphaddot(t)
        alfdd = -amp*(omega^2)*sin(omega*t)
        return alfdd*(pi/180)
    end
    
    return U, Udot, alpha, alphadot, alphaddot
end

nothing

function states!(dx, x, p, t)
    ### Unpack p
    V, Vdot, alpha, alphadot, c, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv = p

    
    ### Environmental values
    M = V(t)/a #Mach number
    beta = sqrt(1-(M^2)) #Compressibility factor #TODO: This will need to be updated to handle if someone hands it a mach number greater than 1. 
    q(t) = alphadot(t)*c/V(t) #Pitch rate

    ### Constants
#     conv = (c/(2*V(t)))
    tau_p = T_p*conv(t)
    tau_f = T_f*conv(t)
    tau_vl = T_vl*conv(t) #Todo: Should I multiply this by conv
    tau_v = T_v*conv(t)

    # println(tau_vl)

    tau_l = c/a #Todo. I need to find out if this should be the definition just after equation 4 in leishman 1990, or on the top right page of 839 of the same paper. This is defined in the nomenclature section. 
    k_alpha = 1/((1-M) + (pi*beta*(M^2)*((A[1]*b[1])+(A[2]*b[2]))))
    k_q = 1/((1-M) + (2*pi*beta*(M^2)*((A[1]*b[1])+(A[2]*b[2]))))

    ### Update the attached flow states (The first four)
    dx[1] = (2*V(t)/c)*(beta^2)*(-b[1]*x[1]) + alpha(t) + q(t)/2 #High frequency shed vortex
    dx[2] = (2*V(t)/c)*(beta^2)*(-b[2]*x[2]) + alpha(t) + q(t)/2 #Low frequency shed vortex
    dx[3] = (-1/(k_alpha*tau_l))*x[3] + alpha(t) #Impulse due to angle of attack
    dx[4] = (-1/(k_q*tau_l))*x[4] + q(t) #Impulse due to pitch rate. 

    ### Calculate the attached lift
    Ccn = dCndalpha*(2*V(t)/c)*(beta^2)*((A[1]*b[1]*x[1])+(A[2]*b[2]*x[2])) #Circulatory lift
    Cina = (4/M)*dx[3] #Normal force due to angle of attack #Todo: Should I calculate this before I update dx? 
    Cinq = (1/M)*dx[4] #Normal force due to pitch rate #Todo: Same as above. I can look at the indicial approach to get an idea. 
    Cpn = Ccn + Cina + Cinq #Attached flow lift
    # println("Cpn: ", Cpn)

    ### Calculate the dynamic flow states (the second four)
    dx[5] = (-1/tau_p)*x[5] + (Cpn/tau_p) #Delayed attached lift.
    # println("tau_p: ", tau_p)

    alpha_f = x[5]/dCndalpha #Equivalent angle of attack
    fp = ffun2(alpha_f, alpha1, S.*conv(t)) #Equivalent separation point #S.*conv(t)
    # println("fp: ", fp)
    # if fp<0
    #     println("")
    #     println("alpha_f: ", alpha_f)
    #     println("fp: ", fp)
    # end
    dx[6] = (-1/tau_f)*x[6] + (fp/tau_f) #Delayed equivalent separation point

    dx[7] = (V(t)/(3*c))*Heavi(x[5]-Cn1)*Heavi(2*tau_vl - x[7]) #- V(t)*Heavi(Cn1-x[5])*Heavi(x[7])/c #Position of the separation point
    if x[5]<Cn1
        x[7]=0
    end
    
    Ccndot = dCndalpha*(2*Vdot(t)/c)*(1-(M^2))*((A[1]*b[1]*x[1])+(A[2]*b[2]*x[2])) + dCndalpha*(2*V(t)/c)*(1-(M^2))*((A[1]*b[1]*dx[1])+(A[2]*b[2]*dx[2])) - dCndalpha*(4*(M^2)/c)*Vdot(t)*((A[1]*b[1]*x[1])+(A[2]*b[2]*x[2])) #Derivative of the circulatory normal force with respect to time. 
    Cvdot = (Ccndot*(1-((1+sqrt(x[6]))^2)/4) - Ccn*(1 + sqrt(x[6]))*dx[6]/(4*sqrt(x[6])))*Heavi(2*tau_vl - x[7])*Heavi(x[5]-Cn1) #Derivative of the vortex lift contribution with respect to time.
    
#     if Cvdot>0.001
#         println(Cvdot)
#     end
    # println("") 
    # println("Ccn: ", Ccn)
    # println("Ccndot: ", Ccndot)
    # println("Cvdot: ", Cvdot)
    dx[8] = (-1/tau_v)*x[8] + (Cvdot/tau_v) #Vortex lift contribution

end

function parsesolution2(sol, p; eta=0.95)
    function extractdata(sol)
        t = sol.t
        u = reduce(hcat, sol.u)'
        n,m = size(u)
        du = zeros(n,m)

        for j=1:m
            du[:,j] = gradient(t, u[:,j], t)
        end
        return t, u, du
    end

    V, Vdot, alpha, alphadot, c, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv = p

    t, x, dx = extractdata(sol)

    n = length(t)

    #Initialize arrays
    Ccn = zeros(n)
    Cina = zeros(n)
    Cinq = zeros(n)
    Cpn = zeros(n)

    alphaE = zeros(n)
    Cnf = zeros(n)
    Cnt = zeros(n)
    Cnt2 = zeros(n)
    
    Cfc = zeros(n)
    Cdn = zeros(n)

    for i = 1:n
        ### Calculate constants
        M = V(t[i])/a
        beta = sqrt(1-(M^2))

        ### Calculate attached lift
        Ccn[i] = dCndalpha*(2*V(t[i])/c)*(beta^2)*((A[1]*b[1]*x[i,1])+(A[2]*b[2]*x[i,2])) #Circulatory lift
        Cina[i] = (4/M)*dx[i,3] #Normal force due to angle of attack  
        Cinq[i] = (1/M)*dx[i,4] #Normal force due to pitch rate 

        Cpn[i] = Ccn[i] + Cina[i] + Cinq[i] #Attached normal force

        alphaE[i] = (beta^2)*(2*V(t[i])/c)*((A[1]*b[1]*x[i,1])+(A[2]*b[2]*x[i,2])) #Effective angle of attack

        Cnf[i] = dCndalpha*(((1+sqrt(x[i,6]))^2)/4)*alphaE[i] + Cina[i] + Cinq[i] #Nonlinear dynamic normal force

        Cnt[i] = Cnf[i] + x[i,8] #Total dynamic normal force
        Cnt2[i] = (((1+sqrt(x[i,6]))^2)/4)*(dCndalpha*alphaE[i] + Cina[i] + Cinq[i]) + x[i,8]
        
        Cfc[i] = eta*dCndalpha*(alphaE[i]^2)*sqrt(x[i,6])
        Cdn[i] = Cnt[i]*sin(alpha(t[i])) - Cfc[i]*cos(alpha(t[i]))
    end
    return t, x, dx, Cnt, Cpn, Cdn, Cfc, Cnt2
end

nothing

polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989staticdata_fig1.csv", ',')


polar[:,1] = polar[:,1].*(pi/180) #Convert to radians

zeroalf, zeroidx = nearestto(polar[:,2], 0.0) #Nearest to zero lift
tenalf, tenidx = nearestto(polar[:,1], 10.0*pi/180) #Ten is typically before stall #*pi/180
linrng = zeroidx:tenidx

mxb = fit(polar[linrng,1], polar[linrng,2], 1)

dCndalpha_exp = mxb.coeffs[2] #least squares fit 

alpha0 = roots(mxb)[1] #least squares fit   # NACA 0012 is a symmetric airfoil... so shouldn't the zero lift be at zero? 

# polarplt = plot(polar[:,1], polar[:,2], leg=false)
# plot!(mxb)
# display(polarplt)

maxcn, maxcnidx = findmax(polar[:,2]) #critical lift coefficient (given by Larsen 2007), the beginning of LE separation. 
# println("Static Max Cl: ", maxcl)
alphav = polar[maxcnidx, 1] 

maxcn = dCndalpha_exp*(alphav-alpha0) #Projecting the inviscid lift to the static stall point
Cn1_exp, Cn1idx = nearestto(polar[:,2], maxcn)
alpha1_exp = polar[Cn1idx, 1]

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

function ffun2(alpha, alpha1, S)
    if alpha<=alpha1
        return 1.0-0.3*exp((alpha-alpha1)/S[1])
    else
        return 0.04 + 0.66*exp((alpha1-alpha)/S[2]) #One paper has a plus here, another has a minus. 
    end
end

nothing

c1 = .415

Ufun1, Udotfun1, alphafun1, alphadotfun1, alphaddotfun1 = prepenvironment(; c=c1, M=0.3, a=343.3, amp=10.0, shift=10.0, k=0.1)

conv1(t) = (c1/(2*Ufun1(t)))

p1 = [Ufun1, Udotfun1, alphafun1, alphadotfun1, c1, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv1]

prob1 = ODEProblem(states!, u0, tspan, p1)
sol1 = solve(prob1;dtmax=0.0001)

t1, u1, du1, Cnt1, Cnp1, Cdn1, Cfc1, Cntt1 = parsesolution2(sol1, p1)


alphavec1 = alphafun1.(t1)
alfavec1 = alphavec1.*(180/pi)
nothing

c2 = .415

Ufun2, Udotfun2, alphafun2, alphadotfun2, alphaddotfun2 = prepenvironment(; c=c2, M=0.3, a=343.3, amp=10.0, shift=10.0, k=0.1)

conv2(t) = (c2/(1*Ufun1(t)))

p2 = [Ufun2, Udotfun2, alphafun2, alphadotfun2, c2, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv2]

prob2 = ODEProblem(states!, u0, tspan, p2)
sol2 = solve(prob2;dtmax=0.0001)

t2, u2, du2, Cnt2, Cnp2, Cdrag2, Cfchord2m, Cntt2 = parsesolution2(sol2, p2)


alphavec2 = alphafun2.(t2)
alfavec2 = alphavec2.*(180/pi)
nothing

c3 = .415

Ufun3, Udotfun3, alphafun3, alphadotfun3, alphaddotfun3 = prepenvironment(; c=c3, M=0.3, a=343.3, amp=10.0, shift=10.0, k=0.1)

conv3(t) = (c3/(1*Ufun3(t)))

p3 = [Ufun3, Udotfun3, alphafun3, alphadotfun3, c3, dCndalpha_exp, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv3]

prob3 = ODEProblem(states!, u0, tspan, p3)
sol3 = solve(prob3;dtmax=0.0001)

t3, u3, du3, Cnt3, Cnp3, Cdn3, Cfc3, Cntt3 = parsesolution2(sol3, p3)


alphavec3 = alphafun3.(t3)
alfavec3 = alphavec3.*(180/pi)
nothing

### Try artificially increasing the slope
c4 = .415

Ufun4, Udotfun4, alphafun4, alphadotfun4, alphaddotfun4 = prepenvironment(; c=c4, M=0.3, a=343.3, amp=10.0, shift=10.0, k=0.1)

conv4(t) = (c4/(1*Ufun4(t)))

p4 = [Ufun4, Udotfun4, alphafun4, alphadotfun4, c4, dCndalpha*1.1, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv4]

prob4 = ODEProblem(states!, u0, tspan, p4)
sol4 = solve(prob4;dtmax=0.0001)

t4, u4, du4, Cnt4, Cnp4, Cdn4, Cfc4, Cntt4 = parsesolution2(sol4, p4)


alphavec4 = alphafun4.(t4)
alfavec4 = alphavec4.*(180/pi)
nothing

S = S.*(c2/v)

include("BeddoesLeishmanss.jl")

v = M*a
Tp = T_p*(c2/v)
Tf = T_f*(c2/v)
Tv = T_v*(c2/v)
Tvl = T_vl*(c2/v)



p5 = [Ufun2, alphafun2, alphadotfun2, c2, dCndalpha, alpha1, alpha0, Cn1, A, b, a, Tp, Tf, Tv, Tvl]

prob5 = ODEProblem(states_liftonly!, u0, tspan, p5)
sol5 = solve(prob5;dtmax=0.0001)

Cdn, Cpn, Cdd, Cdd2, Cdd3, Cdd4, Cdd5, Cdd6, Cdd7, Cdd8, Cdd9, Cdd10, Cdd11, Cdd12, Cdd13, Cdd14, Cdd15, Cdd16, Cdd17, Cdd18, Cdd19, u5, du5, Cc, Cfc, Cfc2, Cfc3, Cfc4, Cfc5, Cfc6, Cfc7, Cfc8, Cfc9, Cfc10, Cfc11, Cfc12, Cfc13, Cfc14, Cfc15, Cfc16, Cfc17, Cfc18, Cfc19 = parsesolution(sol5, p5)

alphavec5 = alphafun2.(sol5.t)
alfavec5 = alphavec5.*(180/pi)

nothing

exppolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/leishman1989statespace/fig9/naca0012/experimental data.csv", ',')

modpolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/leishman1989statespace/fig9/naca0012/Model.csv", ',')

nothing

cycleplt = plot(legend=:topleft, xaxis="Angle of Attack (degrees)", yaxis="Normal Coefficient", title="NACA 0012")
scatter!(exppolar[:,1], exppolar[:,2], lab="Experimental")
scatter!(modpolar[:,1], modpolar[:,2], lab="Model - paper")
# plot!(alfavec1, Cnt1, lab="BL c/2V")
plot!(alfavec2, Cnt2, lab="BL c/V")
# plot!(alfavec3, Cnt3, lab="BL c/v, exp slope")
# plot!(alfavec4, Cnt4, lab="BL c/v, 1.1 slope")
# plot!(alfavec4, Cnp4, lab="Cnp4")
plot!(alfavec5, Cdn, lab="First model")
display(cycleplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/comparison/compareblmodels_naca0012_10.0a_10.0s.png")

statesplt2 = plot(sol2,linewidth=2,xaxis="t",title=["x1" "x2" "x3" "x4" "C\'n" "f\'\'" "tau" "Cvn"], leg=false,layout=(2,4))
plot!(sol5,linewidth=2,xaxis="t",title=["x1" "x2" "x3" "x4" "C\'n" "f\'\'" "tau" "Cvn"], leg=false,layout=(2,4))
display(statesplt2)

println("max(U_2): ", maximum(u2[:,5]))
println("max(U_5): ", maximum(u5[:,5]))

println(Tvl)

println(T_vl*(c2/v))

Ufunsc, Udotfunsc, alphafunsc, alphadotfunsc, alphaddotfunsc = prepenvironment(; c=0.1, M=0.379, a=343.3, amp=8.1, shift=10.3, k=0.075)

nothing

S = S./(c2/v)
c6 = 0.1

convsc(t) = (c6/(1*Ufunsc(t)))

p6 = [Ufunsc, Udotfunsc, alphafunsc, alphadotfunsc, c6, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, convsc]

prob6 = ODEProblem(states!, u0, tspan, p6)
sol6 = solve(prob6;dtmax=0.0001)

t6, u6, du6, Cnt6, Cpn6, Cdn6, Cfc6, Cntt6 = parsesolution2(sol6, p6)


alphavec6 = alphafunsc.(t6)
alfavec6 = alphavec6.*(180/pi)

nothing

S = S.*(c6/v)

include("BeddoesLeishmanss.jl")

v = M*a
Tp = T_p*(c6/v)
Tf = T_f*(c6/v)
Tv = T_v*(c6/v)
Tvl = T_vl*(c6/v)



p7 = [Ufunsc, alphafunsc, alphadotfunsc, c6, dCndalpha, alpha1, alpha0, Cn1, A, b, a, Tp, Tf, Tv, Tvl]

prob7 = ODEProblem(states_liftonly!, u0, tspan, p7)
sol7 = solve(prob7;dtmax=0.0001)

Cdn7, Cpn, Cdrag, Cdd2, Cdd3, Cdd4, Cdd5, Cdd6, Cdd7, Cdd8, Cdd9, Cdd10, Cdd11, Cdd12, Cdd13, Cdd14, Cdd15, Cdd16, Cdd17, Cdd18, Cdd19, u5, du5, Cc, Cfc, Cfc2, Cfc3, Cfc4, Cfc5, Cfc6, Cfc7, Cfc8, Cfc9, Cfc10, Cfc11, Cfc12, Cfc13, Cfc14, Cfc15, Cfc16, Cfc17, Cfc18, Cfc19 = parsesolution(sol7, p7)

alphavec7 = alphafunsc.(sol7.t)
alfavec7 = alphavec7.*(180/pi)

nothing

exppolar2 = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989fig8Cn.csv", ',')

cycleplt = plot(legend=:topleft, xaxis="Angle of Attack (degrees)", yaxis="Normal Coefficient", title="NACA 0012")
scatter!(exppolar2[:,1], exppolar2[:,2], lab="Experimental")
plot!(alfavec7, Cdn7, lab="First Model")
plot!(alfavec6, Cnt6, lab="Second Model")
display(cycleplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/comparison/compareblmodels_naca0012_8.1a_10.3s_091621_paperinputs.png")

dragexppolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/leishman1989statespace/fig9/naca0012/leishman1989ss fig9 drag experiment.csv", ',')
dragmodpolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/leishman1989statespace/fig9/naca0012/leishmanss1989 drag naca0012 fig 9 model.csv", ',')
nothing

dragplt = plot(legend=:topleft, xaxis="Angle of Attack (degrees)", yaxis="Coefficient of Drag")
scatter!(dragexppolar[:,1], dragexppolar[:,2], lab="Experimental")
scatter!(dragmodpolar[:,1], dragmodpolar[:,2], lab="Their Model")
plot!(alfavec2, Cdrag2, lab="My 2nd Model")
plot!(alfavec5, Cdd, lab="My 1st Model")
# plot!(alfavec2, Cfchord2, lab="2nd model chordwise force")
# plot!(alfavec5, Cc, lab="1st Model chordwise force")
display(dragplt)

using Statistics
mach = 0.3
a = 343
S = [2.75, 1.4] 

dtminn = 0.002
alg = RK4()

fvec = [4, 5.33, 6, 8, 10]
kcvec = [0.051, 0.067, 0.076, 0.102, 0.127]

ceevec = zeros(length(fvec))

for i = 1:length(fvec)
    ceevec[i] = kcvec[i]*mach*a/(fvec[i]*pi)
end

cee = round(mean(ceevec); sigdigits=4)

Ufun8, Udotfun8, alphafun8, alphadotfun8, alphaddotfun8 = prepenvironment(; c=cee, M=mach, a=343.3, amp=10.0, shift=10.0, k=0.1)

conv8(t) = (cee/(1*Ufun8(t)))

p8 = [Ufun8, Udotfun8, alphafun8, alphadotfun8, cee, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv8]

prob8 = ODEProblem(states!, u0, tspan, p8)
sol8 = solve(prob8, alg; dtmin=dtminn, force_dtmin=true) #;dtmax=0.0001

t8, u8, du8, Cnt8, Cnp8, Cdrag8, Cfchord8 = parsesolution2(sol8, p8)


alphavec8 = alphafun8.(t8)
alfavec8 = alphavec8.*(180/pi)

v = mach*a
S = S.*(cee/v)
Tp = T_p*(cee/v)
Tf = T_f*(cee/v)
Tv = T_v*(cee/v)
Tvl = T_vl*(cee/v)

p9 = [Ufun8, alphafun8, alphadotfun8, cee, dCndalpha, alpha1, alpha0, Cn1, A, b, a, Tp, Tf, Tv, Tvl]

prob9 = ODEProblem(states_liftonly!, u0, tspan, p9)
sol9 = solve(prob9, alg; dtmin=dtminn, force_dtmin=true) #;dtmax=0.0001

Cdn9, Cpn9, Cdrag9, Cdd2, Cdd3, Cdd4, Cdd5, Cdd6, Cdd7, Cdd8, Cdd9, Cdd10, Cdd11, Cdd12, Cdd13, Cdd14, Cdd15, Cdd16, Cdd17, Cdd18, Cdd19, u5, du5, Cc, Cfc, Cfc2, Cfc3, Cfc4, Cfc5, Cfc6, Cfc7, Cfc8, Cfc9, Cfc10, Cfc11, Cfc12, Cfc13, Cfc14, Cfc15, Cfc16, Cfc17, Cfc18, Cfc19 = parsesolution(sol9, p9)

alphavec9 = alphafun8.(sol9.t)
alfavec9 = alphavec9.*(180/pi)
nothing

cycleplt = plot(legend=:topleft, xaxis="Angle of Attack (degrees)", yaxis="Normal Coefficient", title="NACA 0012 - Mach = $mach")
scatter!(exppolar[:,1], exppolar[:,2], lab="Experimental")
scatter!(modpolar[:,1], modpolar[:,2], lab="Model - paper")
plot!(alfavec9[100:end], Cdn9[100:end], lab="1st Model")
plot!(alfavec8[100:end], Cnt8[100:end], lab="2nd Model")
display(cycleplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/comparison/compareblmodels_naca0012_10.0a_10.0s09162021_differentsolvers.png")

dragplt = plot(legend=:topleft, xaxis="Angle of Attack (degrees)", yaxis="Coefficient of Drag")
scatter!(dragexppolar[:,1], dragexppolar[:,2], lab="Experimental")
scatter!(dragmodpolar[:,1], dragmodpolar[:,2], lab="Their Model")
plot!(alfavec8, Cdrag8, lab="My 2nd Model")
plot!(alfavec9, Cdrag9, lab="My 1st Model")
# plot!(alfavec2, Cfchord2, lab="2nd model chordwise force")
# plot!(alfavec5, Cc, lab="1st Model chordwise force")
display(dragplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/comparison/compareblmodels_naca0012_drag_10.0a_10.0s09162021_differentsolvers.png")

statesplt9 = plot(sol9,linewidth=2,xaxis="t",title=["x1" "x2" "x3" "x4" "C\'n" "f\'\'" "tau" "Cvn"], leg=:topright,layout=(2,4), lab="1st model")
plot!(sol8,linewidth=2,xaxis="t",title=["x1" "x2" "x3" "x4" "C\'n" "f\'\'" "tau" "Cvn"], lab="2nd model",layout=(2,4))
display(statesplt9)

mach = 0.3
kay = 0.10
a = 343.3
S = [3.0, 2.3] #[2.75, 1.4] 
dCndalpha = 0.108*(180/pi) #0.113*(180/(1*pi))
alpha1 = 15.25*(pi/180)
Cn1 = 1.45
T_p = 1.7
T_f = 3.0
T_v = 6.0
T_vl = 7.0
a = 343.0

dtminn = 0.0001
alg = Tsit5()

fvec = [4, 5.33, 6, 8, 10]
kcvec = [0.051, 0.067, 0.076, 0.102, 0.127]

ceevec = zeros(length(fvec))

for i = 1:length(fvec)
    ceevec[i] = kcvec[i]*mach*a/(fvec[i]*pi)
end

cee = round(mean(ceevec); sigdigits=4)

Ufun10, Udotfun10, alphafun10, alphadotfun10, alphaddotfun10 = prepenvironment(; c=cee, M=mach, a=343.3, amp=10.0, shift=5.0, k=kay)

conv10(t) = (cee/(1*Ufun10(t)))

p10 = [Ufun10, Udotfun10, alphafun10, alphadotfun10, cee, dCndalpha, alpha1, Cn1, A, b, S, T_p, T_f, T_vl, T_v, a, conv10]

prob10 = ODEProblem(states!, u0, tspan, p10)
sol10 = solve(prob10, alg; dtmin=dtminn, force_dtmin=true) #;dtmax=0.0001

t10, u10, du10, Cnt10, Cnp10, Cdrag10, Cfchord10, Cntt10 = parsesolution2(sol10, p10)


alphavec10 = alphafun10.(t10)
alfavec10 = alphavec10.*(180/pi)

v = mach*a
S = S.*(cee/v)
Tp = T_p*(cee/v)
Tf = T_f*(cee/v)
Tv = T_v*(cee/v)
Tvl = T_vl*(cee/v)

p11 = [Ufun10, alphafun10, alphadotfun10, cee, dCndalpha, alpha1, alpha0, Cn1, A, b, a, Tp, Tf, Tv, Tvl]

prob11 = ODEProblem(states_liftonly!, u0, tspan, p11)
sol11 = solve(prob11, alg; dtmin=dtminn, force_dtmin=true) #;dtmax=0.0001

Cdn11, Cpn11, Cdrag11, Cdd2, Cdd3, Cdd4, Cdd5, Cdd6, Cdd7, Cdd8, Cdd9, Cdd10, Cdd11, Cdd12, Cdd13, Cdd14, Cdd15, Cdd16, Cdd17, Cdd18, Cdd19, u5, du5, Cc, Cfc, Cfc2, Cfc3, Cfc4, Cfc5, Cfc6, Cfc7, Cfc8, Cfc9, Cfc10, Cfc11, Cfc12, Cfc13, Cfc14, Cfc15, Cfc16, Cfc17, Cfc18, Cfc19 = parsesolution(sol11, p11)

alphavec11 = alphafun10.(sol11.t)
alfavec11 = alphavec11.*(180/pi)
nothing

attachedlift = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/leishman1989statespace/fig8/naca0012/leishman1989 fig 8 lift naca0012.csv", ',')
mal = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/leishman1989statespace/fig8/naca0012/leishman1989 fig 8 naca 0012 lift model.csv", ',')

using FLOWMath

liftfit = Akima(polar[:,1], polar[:,2])

Csn11 = liftfit.(alphavec11)

nothing

cyclepltattached = plot(legend=:topleft, xaxis="Angle of Attack (degrees)", yaxis="Normal Coefficient", title="NACA 0012 - Mach = $mach")
scatter!(attachedlift[:,1], attachedlift[:,2], lab="Experimental")
scatter!(mal[:,1], mal[:,2], lab="Model - paper")
plot!(alfavec11, Cdn11, lab="1st Model")
plot!(alfavec11, Cpn11, lab="1st attached") #For some odd reason, the stall model is not firing. 
plot!(alfavec10, Cnt10, lab="2nd Model")
plot!(alfavec11, Csn11, lab="Static")
display(cyclepltattached)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/comparison/compareblmodels_naca0012_10.0a_10.0s09162021_attachedflow.png")

statesplt10 = plot(sol10,linewidth=2,xaxis="t",title=["x1" "x2" "x3" "x4" "C\'n" "f\'\'" "tau" "Cvn"], leg=false,layout=(2,4), lab="2nd model")
plot!(sol11,linewidth=2,xaxis="t",title=["x1" "x2" "x3" "x4" "C\'n" "f\'\'" "tau" "Cvn"], lab="1st model",layout=(2,4))
display(statesplt10)

Cn1

dCndalpha




