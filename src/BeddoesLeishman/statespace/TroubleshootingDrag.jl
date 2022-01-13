using DifferentialEquations
using Roots
using DelimitedFiles
using Plots
using FLOWMath
using Polynomials

function Heavi(t)
    if t>0
        return 1
    else
        return 0
    end
end




"""
    BeddoesLeishman(U, Udot, V, Vdot, alpha, alphadot, alphaddot, c, dcldalpha, alpha0, maxcl, A, S, w; ffunction=NaN, u0=zeros(8), a=343)

### Inputs:
- U::function - This function should return the freestream velocity (X direction) as a function of time (meters/second)
- Udot::function - derivative of U as a function of time (meters/second^2)
- V::function - vertical freestreem velocity as a function of time. Currently unused. (meters/second)
- Vdot::function - derivative of V as a function of time (meters/second^2)
- alpha::function - angle of attack (in radians) as a function of time (radians)
- alphadot::function - time derivative of alpha (radians/second)
- alphaddot::function - second time derivative of alpha (radians/second^2)
- c::Float64 - chord value (meters)
- dcldalpha::Float64 - slope of the static lift curve (Cl/alpha) (1/radians)
- alpha0::Float64 - zero lift angle of attack (radians)
- maxcl::Float64 - separation lift coefficient (often the max value) (unitless)
- A::Array{Float64, 1} - Constants describing how the airfoil responds to pressure changes and implusive changes. Length 4
- S::Array{Float64, 1} - Constants describing how the airfoil stalls
- w::Array{Float64, 1} - 

### Outputs:

### Notes:
"""
function BeddoesLeishman(U, Udot, V, Vdot, alpha, alphadot, alphaddot, c, dcldalpha, alpha0, alphav, maxcl, A, S, w; ffunction=:Leishman, u0=zeros(8), a=343, Cd0=NaN)

    if length(u0)!=8
        error("Length of inital states must be 8.")
    end

    if ffunction==:Leishman
        function ffun(alphan, alpha1, alpha0) #We could try some other f functions.
            if alphan<=alpha1
                return 1.0 - 0.3*exp((alphan-alpha1)/S[1])
            else
                return 0.04 + 0.66*exp((alpha1-alphan)/S[2])
            end
        end
    else
        ffun = ffunction
    end
    
    function alphaefun(aoa, aoadot, c, U) #This function is different when accounting for turbulence and drafts perpendicular to the airfoil. (The same function (EQ 30), I've just ommitted some terms and simplified for 2D.)
        return aoa + c*aoadot/(2*U)
    end

    function CL0(alpha)
        # alpha0 = 0.01579 #-0.019
        # dcldalpha = 6.7219 #2*pi
        return dcldalpha*(alpha-alpha0)
    end

    function CL0dot(alphadot)
        # dcldalpha = 6.7219 #2*pi
        return dcldalpha*alphadot
    end

    p = [CL0, CL0dot, dcldalpha, maxcl, U, alpha, alphadot, alphaddot, alphav, alpha0, A, S, w, c, a, Cd0]
    #TODO: Maybe I'll make it so the function creates u0 based off of the initial values of the u somehow.... don't know how, but somehow. 
    
    function states!(du, u, p, t)
    
        #Cl0dot - input function - So I think this is the time derivative of the inviscid lift coefficient... Which.... isn't part of the ode.... but changes with time. I should be able to take the derivative. I can, so I can just pass that function in, or whatever. 
        #theta -  eq 7 inverse
        #alphav - separation angle
        #alpha - Input function
        #alphadot - Input function
        CL0, CL0dot, dcldalpha, maxCl, U, alpha, alphadot, alphaddot, alpha1, alpha0, A, S, w, c, a, Cd0 = p
    
        CpL0v = maxCl #Todo: I currently have this as the max Cl, but maybe it should be higher? Or lower? 
        M = U(t)/a
        alphae = alphaefun(alpha(t), alphadot(t), c, U(t))
        # println(alpha(t))
        # println(alphae)
        # println("")
        
        du[1] = -w[1]*u[1] + A[1]*CL0dot(alphae) #This should be found using alphaE #Checked against EQ A.5
        du[2] = -w[2]*u[2] + A[2]*CL0dot(alphae) #This should be found using alphaE #Checked against EQ A.5
        du[3] = -w[5]*u[3] + 4*A[3]*alphadot(t)/M #Checked against EQ A.5
        du[4] = -w[6]*u[4] + A[4]*c*alphaddot(t)/(M*U(t)) #Todo: I should check that the V in EQ A.5 is indeed the freestream velocity. 
        du[5] = w[7]*(CL0(alpha(t)) - u[1] - u[2] + u[3] + u[4] - u[5]) #Note: Switched to using the alpha input function. -> Made no difference. #Checked against EQ A.5
    
        alphaf = u[5]/dcldalpha + alpha0 #Todo. Do we want alpha f or alpha e? -> Should be alpha f, according to EQ A.5 # Checked against the equation given in the paragraph below A.3
        du[6] = w[3]*(ffun(alphaf, alpha1, alpha0) - u[6]) #Note:  #Currently using the f function from the Beddoes-Leishman explanation, because that seemed reasonable. #X
    
        ### u[7] is the position of the LE vortex. thus du[7] is the velocity of the LE vortex. I have conditions to halt the vortex, and zero the vortex. Probably not continuous. 
        if u[7]<1.0
            du[7] = U(t)*Heavi(u[5] - CpL0v)/(3*c)
        else
            du[7] = 0 
        end
        if (u[5] - CpL0v)<0  
            u[7] = 0
        end
    
        CL0d = CL0(alpha(t)) - u[1] - u[2]
        CL0ddot = CL0dot(alphadot(t))-du[1]-du[2] #Note: I changed this to the derivative of Cl0 using alphadot instead of t. It made a small difference. 
        # delCldot = (1 - ((u[6]-0.5)/2)^4)*CL0ddot - 4*(((u[6]-0.5)/2)^3)*CL0d #Todo. This should probably be rechecked. I don't think that I'm doing this quite right. Rather, I might be doing it right, but it seems weird that I should go analytically calculate it. It means that I wouldn't be able to use it when I don't have a sinusoidal input. -- If so, wouldn't my alpha(t) function have a cosine instead of a sine. .... Oh... it does. 
        # Okay, I've redone it, but with the f formulation (EQ 3). 
        delCldot = (-du[6]/(4*sqrt(u[6])) - du[6]/4)*CL0d + (3/4 - sqrt(u[6])/2 - u[6]/4)*CL0ddot #That made a significant difference. Also, changing the initial condition to it's minimum once it hit it's regular state. I think there is still some kind of problem. I mean, it could be a difference in integration scheme, but I doubt it. Last time I played with that, it didn't really make much of a difference.
        # println(du[6])
        # println(delCldot)
    
        du[8] = -w[4]*u[8] + delCldot*Heavi(1-u[7])*Heavi(alphadot(t)) #Checked against EQ A.5
    end
    return states!, p, u0
end

function parsesolution(sol, p)
    u = reduce(hcat, sol.u)'
    n,m = size(u)
    
    CL0, CL0dot, dcldalpha, maxCl, U, alpha, alphadot, alphaddot, alpha1, alpha0, A, S, w, c, a, Cd0 = p

    function alphaefun(aoa, aoadot, c, U) #This function is different when accounting for turbulence and drafts perpendicular to the airfoil. (The same function (EQ 30), I've just ommitted some terms and simplified for 2D.)
        return aoa + c*aoadot/(2*U)
    end
    
    Cl = zeros(n)
    Cc = zeros(n)
    Cd = zeros(n)
    
    for i = 1:n
        t = sol.t[i]
        
        Cl[i] = (((1 + sqrt(u[i,6]))/2)^2)*(CL0(alpha(t)) - u[i,1] - u[i,2] + u[i,3]+ u[i,4]) + u[i,8] #Todo. Need to double check this. done. Matches Equation A.7

        M = U(t)/a
        beta = sqrt(1-(M^2))
        alphaequiv = (beta^2)*(A[1]*w[1]*u[i,1] + A[2]*w[2]*u[i,2]) 
        # println("alphaequiv 1: ", alphaequiv)
        #Leishmean 1990 EQ 22
        alphaequiv = alphaefun(alpha(t), alphadot(t), c, U(t))
        # println("alphae 2: ", alphaequiv)
        Cc[i] = dcldalpha*(alphaequiv^2) #Leishman 1990 EQ 23
        Cd[i] = Cl[i]*cos(alpha(t))*sin(alpha(t)) - Cc[i]*cos(alpha(t)) #Cl[i]*sin(alpha(t)) + Cc[i]*cos(alpha(t)) #Leishman 1990 EQ 24 
        #Todo. Shouldn't this have the zero lif drag? -> I think so, because this should be the pressure drag, and the zero lift drag would be the friction drag. 
        #Least of my worries apparently, It appears to need something else. Or maybe I misunderstand something. 
        if !isnan(Cd0)
            Cd[i] += Cd0
            println("Added friction drag")
        end
        # println("")
        # println("M: ", M)
        # println("beta: ", beta)
        # println("alphaequiv: ", alphaequiv)
        # println("Cc: ", Cc[i])
        # println("Cd: ", Cd[i])
        # println("")
        # println("")
        
    end
    return Cl, Cd, Cc
end

function U(t)
    M = 0.379
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
    M = 0.379
    a = 343.3
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
    a = 343.3
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
    a = 343.3
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
M = 0.379
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
a = 343.0 #Todo: I'm not entirely sure what this will do to the computations, but I think it'll have an effect. I'm not sure what to set it. 

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
nothing

function alphaefun(aoa, aoadot, c, U) #This function is different when accounting for turbulence and drafts perpendicular to the airfoil. (The same function (EQ 30), I've just ommitted some terms and simplified for 2D.)
    return aoa + c*aoadot/(2*U)
end

alphavec = alpha.(sol.t)
u = reduce(hcat, sol.u)'
n,m = size(u)
    
CL0, CL0dot, dcldalpha, maxCl, U, alpha, alphadot, alphaddot, alpha1, alpha0, A, S, w, c, a, Cd0 = p

    
    
Cl = zeros(n)
Cc = zeros(n)
Cd = zeros(n)
    
for i = 1:n
    t = sol.t[i]
        
    Cl[i] = (((1 + sqrt(u[i,6]))/2)^2)*(CL0(alpha(t)) - u[i,1] - u[i,2] + u[i,3]+ u[i,4]) + u[i,8] #Todo. Need to double check this. done. Matches Equation A.7
end       






### Here, we're going to test if the output of Larsen's algorithm is indeed lift, or if it is normal force. I should be able to multiply each time step by cos(α) to get the normal coefficient... I think.... At least... that's what the image up above seems to insinuate. (I might have to multiply by cosine alphae)

Cn = zeros(n)
for i = 1:n
    Cn[i] = Cl[i]*cos(alphavec[i])
end


exppolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989fig8Cn.csv", ',')

inviscid(alf) = dcldalpha.*(alf.-alpha0)
cli = inviscid(alphavec)

cycleplot = plot(alpha.(sol.t).*(180/pi), Cl, lab="BLL", leg=:topleft)
plot!(alpha.(sol.t).*(180/pi), Cn, lab="Cn")
# plot!(alpha.(sol.t).*(180/pi), cli, lab="inviscid")
# plot!(aoavec.*(180/pi), clo, lab="inviscid")
scatter!(exppolar[:,1], exppolar[:,2], lab="experimental")
xlabel!("angle (degrees)")
ylabel!("Coefficient of Lift")
# ylims!((0.0,1.5))
# display(cycleplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanlarsen/leishman1989fig8/NACA0012_forwardfig8_experimentalslope_mintimestep0.0001.png")



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

### Geometry
c = 0.1
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

alphavec = alpha.(sol.t)
u = reduce(hcat, sol.u)'
n,m = size(u)
    
CL0, CL0dot, dcldalpha, maxCl, U, alpha, alphadot, alphaddot, alpha1, alpha0, A, S, w, c, a, Cd0 = p

Cl = zeros(n)
    
for i = 1:n
    t = sol.t[i]
        
    Cl[i] = (((1 + sqrt(u[i,6]))/2)^2)*(CL0(alpha(t)) - u[i,1] - u[i,2] + u[i,3]+ u[i,4]) + u[i,8] #Todo. Need to double check this. done. Matches Equation A.7
end  



### Here, we're going to test if the output of Larsen's algorithm is indeed lift, or if it is normal force. I should be able to multiply each time step by cos(α) to get the normal coefficient... I think.... At least... that's what the image up above seems to insinuate. (I might have to multiply by cosine alphae)

Cn = zeros(n)
for i = 1:n
    Cn[i] = Cl[i]*cos(alphavec[i])
end


exppolar2 = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/fig5rotated.csv", ',')

cycleplot = plot(alpha.(sol.t).*(180/pi), Cl, lab="BLL", leg=:topleft)
plot!(alpha.(sol.t).*(180/pi), Cn, lab="Cn")
# plot!(alpha.(sol.t).*(180/pi), cli, lab="inviscid")
# plot!(aoavec.*(180/pi), clo, lab="inviscid")
scatter!(exppolar2[:,1], exppolar2[:,2], lab="experimental")
xlabel!("angle (degrees)")
ylabel!("Coefficient of Lift")
# ylims!((0.0,1.5))
# display(cycleplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanlarsen/leishman1989fig8/NACA0012_forwardfig8_experimentalslope_mintimestep0.0001.png")

Cc1 = zeros(n)
Cc2 = zeros(n)
Cd1a = zeros(n)
Cd1b = zeros(n)
Cd1c = zeros(n)
Cd1d = zeros(n)
Cd1e = zeros(n)

Cd2a = zeros(n)
Cd2b = zeros(n)
Cd2c = zeros(n)
Cd2d = zeros(n)
Cd2e = zeros(n)

for i = 1:n 
    t = sol.t[i]
    M = U(t)/a
    beta = sqrt(1-(M^2))
    alphaequiv1 = (beta^2)*(A[1]*w[1]*u[i,1] + A[2]*w[2]*u[i,2]) 
    # println("alphaequiv 1: ", alphaequiv)
    #Leishmean 1990 EQ 22
    
    
    Cc1[i] = dcldalpha*(alphaequiv1^2) #Leishman 1990 EQ 23
    
    Cd1a[i] = Cn[i]*sin(alpha(t)) - Cc1[i]*cos(alpha(t)) #Leishman 1990 EQ 24 
    Cd1b[i] = Cl[i]*sin(alpha(t)) - Cc1[i]*cos(alpha(t)) #Assuming Larsen actually gave us the normal force
    Cd1c[i] = Cl[i]*sin(alpha(t)) #Assuming Larsen actually gave us the normal force - not accounting the chordwise force
    Cd1d[i] = Cl[i]*cos(alpha(t))*sin(alpha(t)) - Cc1[i]*cos(alpha(t))  #This is assuming that alpha_e = alpha - should be the same as curve 1a. 
    Cd1e[i] = Cl[i]*cos(alphaequiv1)*sin(alpha(t)) - Cc1[i]*cos(alpha(t))
    
     
    
    alphaequiv2 = alphaefun(alpha(t), alphadot(t), c, U(t))
    Cc2[i] = dcldalpha*(alphaequiv2^2) #Leishman 1990 EQ 23
    Cd2a[i] = Cn[i]*sin(alpha(t)) - Cc2[i]*cos(alpha(t)) #Leishman 1990 EQ 24 
    Cd2b[i] = Cl[i]*sin(alpha(t)) - Cc2[i]*cos(alpha(t))
    Cd2c[i] = Cl[i]*sin(alpha(t))
    Cd2d[i] = Cl[i]*cos(alpha(t))*sin(alpha(t)) - Cc2[i]*cos(alpha(t))
    Cd2e[i] = Cl[i]*cos(alphaequiv2)*sin(alpha(t)) - Cc2[i]*cos(alpha(t))

    println("")
    println(Cd1a[i])
    println(Cd1d[i])
end

alfvec = alpha.(sol.t).*(180/pi)
plt1 = scatter(alfvec, Cd1a, legend=:bottomright, label="Original - Normal - a")
plot!(alfvec, Cd1b, label="Original - Lift - b")
plot!(alfvec, Cd1c, label="Lift*sine(alpha) - c")
xlabel!("angle (degrees)")
ylabel!("Coefficient of Pressure Drag")
display(plt1)


# Cdplot = plot(alpha.(sol.t).*(180/pi), Cdd, legend=:bottomright, label="cosine*sine of lift")
# xlabel!("angle (degrees)")
# ylabel!("Coefficient of Pressure Drag")
# display(Cdplot)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanlarsen/leishman1989fig8/NACA0012DRAG_forwardfig8_experimentalslope_mintimestep0.0001.png")

nothing


plt2 = scatter(alfvec, Cd2a, legend=:topleft, label="Original - Normal - a")
plot!(alfvec, Cd2b, label="Original - Lift - b")
plot!(alfvec, Cd1c, label="Lift*sine(alpha) - c")
xlabel!("angle (degrees)")
ylabel!("Coefficient of Pressure Drag")
display(plt2)

plt3 = scatter(alfvec, Cd1a, legend=:bottomleft, label="Original - Normal - a")
plot!(alfvec, Cd1d, label="Should be same as a - d")
plot!(alfvec, Cd1e, label="Original - Lift - e")
xlabel!("angle (degrees)")
ylabel!("Coefficient of Pressure Drag")
display(plt3)

plt4 = scatter(alfvec, Cd2a, legend=:bottomleft, label="Original - Normal - a")
plot!(alfvec, Cd2d, label="Should be same as a - d")
plot!(alfvec, Cd2e, label="Original - Lift - e")
xlabel!("angle (degrees)")
ylabel!("Coefficient of Pressure Drag")
display(plt4)


