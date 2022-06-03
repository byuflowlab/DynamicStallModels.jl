
function seperationpoint(alpha, afm, afp, clfit, dcldalpha, alpha0)
    #TODO: I'm not really sure that using the minimum of these two is really the way to avoid the problem of this blowing up to infinity. (When alpha=alpha0) (This check happens in the if statement.)

    #TODO: I'm not sure that using the absolute value here is the correct way to dodge the problem of crossing the x axis at different times.

    #Todo. I don't like using if statements for angles of attack outside of the range of attached flow. I want the function to just drive to zero. -> Todo. Plot the seperation function as a function of alpha. -> Removed the if statements and the function naturally drives to zero like the paper (not at the angles as described by Hansen.)
    
    
    ### The theory says that if the angle of attack is less (or greater than) than the seperation point aoa, then the function should return 0. -> However, using these if statements create discontinuities in the seperation point function that don't appear to be in Hansen's implementation. -> It only creates discontinuties if the incorrect values for afm and afp are used. If afm and afp are really where f(alpha)=0 then, the function won't have any discontinuities. 
    # if alpha<afm # -> TODO: If I'm getting rid of this, then I can get rid of the afm and afp arguments. 
    #     return typeof(alpha)(0)
    # elseif alpha>afp
    #     return typeof(alpha)(0)
    # end

    if !(afm < alpha < afp) #Check if alpha is in the bounds. 
        # println("f was set to zero. ")
        return typeof(alpha)(0)
    end
    
    cl_static = clfit(alpha)
    cl_linear = dcldalpha*(alpha-alpha0)
    f = (2*sqrt(abs(cl_static/cl_linear))-1)^2

    if f>1 #Question: What if I don't return this? I might get Inf.... or possibly NaN... but I will less likely get 1.0... which is my problem child in the seperated coefficient of lift function. -> I fixed the fully seperated coefficient of lift function... I just plugged this function inside the other and simplified. 
        return typeof(alpha)(1)
    elseif isnan(f)
        # println("f return NaN")
        return typeof(alpha)(1)
    end

    #Todo. Hansen must have some sort of switch that stops this function from reattaching when the aoa gets really high. -> like the one where you automatically set f=0 when you're outside the bounds of afm, afp
    return f
end



function riso_state_rates!(dx, x, U, Udot, alpha, alphadot, c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp)

    # @show alpha
    Tu = c/(2*U)
    alpha_E = alpha*(1 - A1 - A2) + x[1] + x[2]
    alpha_f = x[3]/dcldalpha + alpha0
    # @show alpha_f
    fst = seperationpoint(alpha_f, afm, afp, clfit, dcldalpha, alpha0)

    dx[1] = -x[1]*(b1 + c*Udot/(2*(U^2)))/Tu + b1*A1*alpha/Tu
    dx[2] = -x[2]*(b2 + c*Udot/(2*(U^2)))/Tu + b2*A2*alpha/Tu
    dx[3] = -x[3]/Tp + (dcldalpha*(alpha_E - alpha0) + pi*Tu*alphadot)/Tp
    dx[4] = -x[4]/Tf + fst/Tf
end

function riso_ode!(dx, x, p, t)
    U, Udot, alpha, alphadot, c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp = p

    riso_state_rates!(dx, x, U(t), Udot(t), alpha(t), alphadot(t), c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp)
end

function static_riso_dae!(resids, dx, x, p, t)
    U, Udot, alpha, alphadot, c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp = p
    @show typeof(x)

    riso_state_rates!(dx, x, U(t), Udot(t), alpha(t), alphadot(t), c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp)

    resids[1] = dx[1]
    resids[2] = dx[1]
    resids[3] = dx[3]
    resids[4] = dx[4]
end

function Clfs(alpha, clfit, dcldalpha, alpha0, afm, afp)
    fst = seperationpoint(alpha, afm, afp, clfit, dcldalpha, alpha0)
    
    clst = clfit(alpha)
    # @show fst
    # if fst>=1 #Todo. This is supposed to happen automagically. 
    #     return clst/2
    # else
    #     return (clst - dcldalpha*(alpha - alpha0)*fst)/(1-fst)
    # end  
    return (clst - dcldalpha*(alpha - alpha0)*fst)/(1-fst)
end

"""
    Clfs(alpha, clfit, dcldalpha, alpha0)

The fully seperated coefficient of lift. The function has the seperation point function inserted, and simplified. 

## Inputs
- alpha - Current angle of attack (radians)
- clfit - A spline object to return the static coefficient of lift as a function of alpha
- dcldalpha - The slope of the attached region of the coefficient of lift vs alpha curve (1/radians). 
- alpha0 - The zero lift angle of attack (radians)
"""
function Clfs(alpha, clfit, dcldalpha, alpha0)
    cl = clfit(alpha)
    a = alpha - alpha0
    am = a*dcldalpha
    term = am*sqrt(abs(cl/(am))) - 3*cl
    return -am*term/(4*cl)
end

function riso_coefficients(x, U, alpha, alphadot, c, clfit, cdfit, dcldalpha, Cd0, alpha0, A1, A2, afm, afp)
    #Todo: Consider creating functions to get the lift, drag and moment seperately. -> Although, the drag depends on the lift. Same with the moment. 

    Tu = c/(2*U)
    alpha_E = alpha*(1 - A1 - A2) + x[1] + x[2]

    # alpha_f = x[3]/dcldalpha + alpha0
    # fst = seperationpoint(alpha_f, afm, afp, clfit, dcldalpha, alpha0)

    ### Calculate Clfs #Todo: Should Clfs use fst calculated from alpha_f or alpha? I'm going to assume it is based off of alpha, but recognize that this could be wrong. 
    # cl_fs = Clfs(alpha_E, clfit, dcldalpha, alpha0, afm, afp)
    cl_fs = Clfs(alpha_E, clfit, dcldalpha, alpha0)
    # @show cl_fs

    Cl_dyn = dcldalpha*(alpha_E - alpha0)*x[4] + cl_fs*(1-x[4]) + pi*Tu*alphadot

    cdst = cdfit(alpha_E)
    fst = seperationpoint(alpha_E, afm, afp, clfit, dcldalpha, alpha0)
    fpp = x[4] #The delayed (unsteady) seperation point. It should never be zero, but perhaps the solver will attempt to use that value. 
    if fpp<0 #Todo: Will this introduce any discontinuities? 
        # println("set fourth state to zero. ") #Note: This is occuring in the converging to static solve. 
        fpp=0
    end
    t1 = (sqrt(fst) - sqrt(fpp))/2 #TODO: Here is the square root that causes problems when the fourth state is negative. 
    t2 = (fst - x[4])/4
    Cd_dyn = cdst + (alpha - alpha_E)*Cl_dyn + (cdst - Cd0)*(t1 - t2)
    return Cl_dyn, Cd_dyn
end

function parsesolution(sol, p, cdfit, Cd0)

    ### Unpack
    x = Array(sol)'
    t = sol.t
    U, Udot, alpha, alphadot, c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp = p

    nt = length(t)

    clvec = zeros(nt)
    cdvec = zeros(nt)

    ### run through the time steps and calculate the dynamic lift and drag (based on the states)
    for i = 1:nt
        ti = t[i]
        clvec[i], cdvec[i] = riso_coefficients(x[i,:], U(ti), alpha(ti), alphadot(ti), c, clfit, cdfit, dcldalpha, Cd0, alpha0, A1, A2, afm, afp)
    end

    return clvec, cdvec, t
end

function parsesolution(sol, p, polar)
    u = reduce(hcat, sol.u)'
    n,m = size(u)
    
    c, A, b, Tp, Tf, dcldalpha, liftfit, U, Udot, V, alpha, alphadot, alpha0 = p

    dragfit = Akima(polar[:,1], polar[:,3])
    momentfit = Akima(polar[:,1], polar[:,4])

    ### Find a^st, the distance between the center of pressure and the quarter chord. #Todo: Use this to get coefficient of moment in riso_coefficients. 
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

function find_seperation_alpha(liftfit, dcldalpha, alpha0; n=10)

    ### Create a residual function to solve. 
    residual(alpha) = abs(dcldalpha*(alpha-alpha0)/4) - abs(liftfit(alpha))

    ### Find a bracketing range for the positive fully stalled angle. 
    aoa = range(0,pi/2; length = n)
    rng = zeros(2)
    for i = 1:n-1
        if residual(aoa[i])*residual(aoa[i+1])<0
            rng[:] = aoa[i:i+1]
            break
        end
    end
    
    ### Find the positive fully stalled angle
    alpha_positive = find_zero(residual, rng, Bisection())

    ### Bracket the negative fully stalled. 
    aoa = range(0,-pi/2; length = n)
    for i = 1:n-1
        if residual(aoa[i])*residual(aoa[i+1])<0
            rng[:] = aoa[i:i+1]
            break
        end
    end

    ### Find the negative fully stalled angle
    alpha_negative = find_zero(residual, rng, Bisection())
    return alpha_positive, alpha_negative
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