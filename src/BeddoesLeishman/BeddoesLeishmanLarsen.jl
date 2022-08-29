using DifferentialEquations
using Roots
using FiniteDiff

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