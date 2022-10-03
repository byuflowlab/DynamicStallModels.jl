# using FLOWMath

export Oye

"""
    Oye(detype::DEType, n::Int, airfoils::Array{Airfoil, 1})

The Øye model struct. It stores airfoil data for every section to be simulated. It can be used as a method to return updated states or state rates depending on it's DEType. 

### Inputs
- detype - The type of model it is, Functional(), Iterative(), or Indicial().
- n - The number of sections to be simulated. 
- airfoils - A vector of Airfoil structs, one corresponding to each section to be simulated. 
"""
struct Oye{TI} <: DSModel
    detype::DEType 
    n::TI #Number of airfoils simulated
    airfoils::Array{Airfoil,1}
end

function nearestto(x, xp)
    residuals = abs.(x .-xp)
    min, idx = findmin(residuals)
    return x[idx], idx
end

function fst(alpha)
    return (2*sqrt(cl(alpha)/clinv(alpha)) - 1)^2
end

function clfs(alpha)
    return (cl(alpha)-(clinv(alpha)*fst(alpha)))/(1-fst(alpha))
end

function tau(c, Vrel; A=4) 
    return A*c/Vrel
end

function fs1(fs, alpha, c, Vrel; delt=0.1)
    return min(fst(alpha) + (fs - fst(alpha))*exp(-delt/tau(c, Vrel)),1.0)
end

function Cl(fs, alpha, c, Vrel)
    fs1_val = fs1(fs, alpha, c, Vrel)
    return fs1_val*clinv(alpha) + (1-fs1_val)*clfs(alpha)
end

function Oye(tspan, alpha, U, c; nt=100, f0=NaN)

    tvec = collect(range(tspan[1], tspan[2]; length=nt))
    fvec = zeros(nt)
    Clvec = zeros(nt)

    deltat = tvec[2]-tvec[1]

    if isnan(f0)
        fvec[1] = fst(alpha(tvec[1]))
    else
        fvec[1] = f0
    end
    # println(fvec[1])
    
    Clvec[1] = Cl(fvec[1], alpha(tvec[1]), c, U(tvec[1]))
    for i = 2:nt
        fvec[i] = fs1(fvec[i-1], alpha(tvec[i]), c, U(tvec[i]); delt=deltat)
        Clvec[i] = fvec[i]*clinv(alpha(tvec[i])) + (1-fvec[i])*clfs(alpha(tvec[i]))
    end

    return Clvec, tvec, fvec
end