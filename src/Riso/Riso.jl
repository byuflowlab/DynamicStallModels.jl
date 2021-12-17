using DifferentialEquations
using FLOWMath
using Roots

"""
6/14/2021 Adam Cardoza
Trying to implement the Riso (Hansen 2004) Dynamic stall model. 

11/1/21 Adam Cardoza
Revisiting to use this instead of the Beddoes-Leishman, because it should be easier. 
"""

# p = [c, A, b, Tp, Tf, C_lalpha, liftfit, U, Udot, alpha, alphadot, alpha0]


function fst(alpha, liftfit, dcdalpha, alpha0)
    # println(alpha)
    # println(alpha0)
    # println(dcdalpha)
    f =  (2*sqrt(abs(liftfit(alpha)/(dcdalpha*(alpha-alpha0)))) - 1)^2
    if f>= 1 || isnan(f)
        return 1.0
    else
        return f
    end
    
    #Todo: I'm not really sure that using the minimum of these two is really the way to avoid the problem of this blowing up to infinity. (When alpha=alpha0)
    #Todo: I'm not sure that using the absolute value here is the correct way to dodge the problem of crossing the x axis at different times. 
end

function Clfs(alpha, liftfit, dcldalpha, alpha0)
    f = fst(alpha, liftfit, dcldalpha, alpha0)
    Cl = 0.0
    if f>= 1.0 #isapprox(f, 1.0, atol=1e-2) ||
        Cl =  liftfit(alpha)/2
    else
        Cl = (liftfit(alpha) - (dcldalpha*(alpha-alpha0))*f)/(1-f) #f cannot be equal to 1. Which it is regularly equal to 1. 
    end
    if isnan(Cl) #TODO: I don't think that I need this anymore. 
        Cl = liftfit(alpha)/2
    end
    # println(Cl)
    return Cl
end

function states!(dx, x, p, t)
    c, A, b, Tp, Tf, dcldalpha, liftfit, U, Udot, V, alpha, alphadot, alpha0 = p
    
    # vbar = V(t) + alphadot(t)*0.5*c #??? Is this how I'd get the downwash at the 3/4 chord? 
    # a34 = vbar/U(t)
    a34 = alpha(t)
    ae = a34*(1-A[1]-A[2]) + x[1] + x[2] #TODO: Could this possibly, somehow be alpha acting on (1-A[1] - A[2]) ? Like, I don't think so, but I guess it always could be. 
    Tu = c/(2*U(t))
    
    dx[1] = (b[1]*A[1]*a34/Tu) - x[1]*(b[1]+ (c*Udot(t)/(2*(U(t))^2)))/Tu 
    dx[2] = (b[2]*A[2]*a34/Tu) - x[2]*(b[2]+ (c*Udot(t)/(2*(U(t))^2)))/Tu
    dx[3] = (dcldalpha*(ae-alpha0) + pi*Tu*alphadot(t))/Tp - x[3]/Tp
    alphaf = (x[3]/dcldalpha)+alpha0
    fp = fst(alphaf, liftfit, dcldalpha, alpha0)
    # println(fp)
    dx[4] = fp/Tf - x[4]/Tf 
    # println(x)
    # println(dx)
end

function parsesolution(sol, p, polar)
    u = reduce(hcat, sol.u)'
    n,m = size(u)
    
    c, A, b, Tp, Tf, dcldalpha, liftfit, U, Udot, V, alpha, alphadot, alpha0 = p

    dragfit = Akima(polar[:,1], polar[:,3])
    momentfit = Akima(polar[:,1], polar[:,4])

    ### Find a^st, the distance between the center of pressure and the quarter chord. 
    alphavec = polar[:,1]
    nn = length(alphavec)
    fvec = zeros(nn)
    astvec = zeros(nn)
    for i=1:nn
        fvec[i] = fst(alphavec[i], liftfit, dcldalpha, alpha0)
        astvec[i] = (momentfit(alphavec[i])-momentfit(alpha0))/liftfit(alphavec[i])
    end

    
    mat = hcat(reverse(fvec), reverse(astvec))
    mat = uniquemat!(mat)
    
    affit = Akima(mat[:,1], mat[:,2])


    
    Cl = zeros(n)
    Cd = zeros(n)
    Cm = zeros(n)
    
    for i = 1:n
        t = sol.t[i]

        a34 = alpha(t)
        ae = a34*(1-A[1]-A[2]) + u[i,1] + u[i,2]
        Tu = c/(2*U(t))
        
        clfs = Clfs(ae, liftfit, dcldalpha, alpha0) 
        
        Cl[i] = dcldalpha*(ae-alpha0)*u[i,4] + clfs*(1-u[i,4]) + pi*Tu*alphadot(t)
        fae = fst(ae, liftfit, dcldalpha, alpha0)
        fterm = (sqrt(fae)-sqrt(u[i,4]))/2 - (fae-u[i,4])/4
        Cd[i] = dragfit(ae) + (alpha(t)-ae)*Cl[i] + (dragfit(ae)-dragfit(alpha0))*fterm
        aterm = affit(u[i,4]) - affit(fst(ae, liftfit, dcldalpha, alpha0))
        # println(aterm)
        Cm[i] = momentfit(ae) + Cl[i]*(aterm) - pi*Tu*alphadot(t)/2
    end
    t = sol.t
    return Cl, Cd, Cm, u, t
end

function findsepalpha(liftfit, dcldalpha, alpha0)
    function residual(alpha)
        return abs(dcldalpha*(alpha-alpha0)/4) - abs(liftfit(alpha))
    end
    rng = (0, 30)
    alphas = zeros(2)
    alphas[1] = find_zero(residual, rng, Bisection())
    # println("Got here")
    rng = (-30, 0)
    alphas[2] = find_zero(residual, rng, Bisection())
    return alphas
end

function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

function uniquemat!(mat;column=1)
    n, m = size(mat)
    if column>m
        error("uniquemat!: You must choose a column within the matrix.")
    end
    listofindices = Int[]
    for i=1:length(mat[:,1])
        value = mat[i,column]
        index = findfirst(x->x==value,mat[:,column])
        push!(listofindices, index)
    end
    unique!(listofindices)
    return mat[listofindices,:]
end

