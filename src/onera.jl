export Onera

"""
    Onera(detype::DEType, n::Int, airfoils::Array{Airfoil, 1})

The Onera model struct. It stores airfoil data for every section to be simulated. It can be used as a method to return updated states or state rates depending on it's DEType. 

### Inputs
- detype - The type of model it is, Continuous() or Discrete().
- n - The number of sections to be simulated. 
- airfoils - A vector of Airfoil structs, one corresponding to each section to be simulated. 
"""
struct Onera{TI} <: DSModel
    detype::DEType 
    n::TI #Number of airfoils simulated
    airfoils::Array{Airfoil,1}
end

#Todo: This needs work. 


function states!(dx, x, y, p, t)
    alpha, alphadot, liftfit, liftfitderiv, dcldalpha, A, omega, zeta = p #Todo: Separate out the parameters from the environmental variables. 

    dx[1] = -omega[1]*x[1] + omega[1]*dcldalpha*(alpha(t)) + (1-A[1])*dcldalpha*alphadot(t)
    dx[2] = x[3]

    delCl = dcldalpha*alpha(t) - liftfit(alpha(t))
    delCldot = dcldalpha*alphadot(t) - liftfitderiv(alpha(t))*alphadot(t)

    dx[3] = -2*zeta*omega[2]*x[3] - (omega[2]^2)*(1 + (zeta^2))*x[2] - (omega[2]^2)*(1 + (zeta^2))*(delCl + A[2]*delCldot)
end

function parsesolution(sol, p)
    u = reduce(hcat, sol.u)'
    n,m = size(u)

    alpha, alphadot, liftfit, liftfitderiv, dcldalpha, A, omega, zeta = p #Todo: 

    t = sol.t

    Cl = zeros(n)

    for i=1:n
        Cl[i] = u[i,1] + u[i,2]
    end
    return Cl, t, u
end