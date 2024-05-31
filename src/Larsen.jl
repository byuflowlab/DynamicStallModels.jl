# using DifferentialEquations
# using Roots

export Larsen

"""
    Larsen(detype::DEType, n::Int, airfoils::Array{Airfoil, 1})

The Larsen model struct. It stores airfoil data for every section to be simulated. It can be used as a method to return updated states or state rates depending on it's DEType. 

### Inputs
- detype - The type of model it is, Functional(), Iterative(), or Indicial().
- n - The number of sections to be simulated. 
- airfoils - A vector of Airfoil structs, one corresponding to each section to be simulated. 
"""
struct Larsen{TI} <: DSModel
    detype::DEType 
    n::TI #Number of airfoils simulated
    airfoils::Array{Airfoil,1}
end

function thetafun(t)
    #Solve equation 7.
    aoa = alpha(t) 
    cl = liftfit(aoa)
    cl0 = CL0(t)
    term = cl/cl0
    if term>1
        #Todo. I'm fudging this here. and quite often.  - see below
        # println("Warning: Cl/Cl0 greater than 1: ", term)
        term = 1.0
    end
    return 4*acos((term)^(1/4)) 
end

function thetafun2(t) #Using the residual form of the function changed nothing.... oddly enough. 
    aoa = alpha(t)
    cl = liftfit(aoa)
    cl0 = CL0(t)
    residual(p) = (cos(p/4)^4)*cl0 - cl

    rng = (0, pi)
    theta = find_zero(residual, rng, Bisection())
    return theta
end

function Heavi(t)
    if t>0
        return 1
    else
        return 0
    end
end

function CL0(t)
    return dcldalpha*(alpha(t)-alpha0)
end

function CL0dot(t)
    return dcldalpha*alphadot(t)
end

function states!(du, u, p, t)

    #Cl0dot - input function - So I think this is the time derivative of the inviscid lift coefficient... Which.... isn't part of the ode.... but changes with time. I should be able to take the derivative. I can, so I can just pass that function in, or whatever. 
    #theta -  eq 7 inverse
    #alphav - separation angle
    #alpha - Input function
    #alphadot - Input function
    CL0, CL0dot, alpha, alphadot, alphav, w = p #Todo: I think I changed p to be just environmental variables. So CL0 and CL0dot would.... oh, they're functions, that should really just be calculated as part of what's going on. So I think this whole model needs a rework. 
    
    du[1] = -w[1]*u[1] + A[1]*CL0dot(t) #C1
    du[2] = -w[2]*u[2] + A[2]*CL0dot(t) #C2
    du[3] = -w[3]*u[3] + w[3]*thetafun(t) #theta_d

    CL0d = CL0(t) - u[1] - u[2] 
    CL0ddot = CL0dot(t)-du[1]-du[2]
    delCldot = (1 - (cos(u[3]/4)^4))*CL0ddot + CL0d*du[3]*sin(u[3]/4)*(cos(u[3]/4)^3) #This is correct, I've double checked it. 
    

    du[4] = -w[4]*u[4] + delCldot*Heavi(alphav-alpha(t))*Heavi(alphadot(t)) #Clv
end

function parsesolution(sol, p)
    u = reduce(hcat, sol.u)'
    n,m = size(u)
    
    CL0, CL0dot, alpha, alphadot, alphav, w = p
    
    Cl = zeros(n)
    
    for i = 1:n
        t = sol.t[i]
        
        Cl[i] = (cos(u[i,3]/4)^4)*(CL0(t) - u[i,1] - u[i,2]) + u[i,4]
        
    end
    return Cl
end