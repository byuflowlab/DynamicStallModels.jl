#=


    #X[1] - 
    #X[2] - 
    #X[3] - 
    #X[4] - 
    # u - 
    # theta - Angle of attack

    # dcldalpha = airfoil.dcldalpha
    # alpha0 = airfoil.alpha0
    # liftfit = airfoil.cl
    # A1 = airfoil.A[1]
    # A2 = airfoil.A[2]
    # b1 = airfoil.b[1]
    # b2 = airfoil.b[2]
    # Tp = airfoil.T[1]
    # Tf = airfoil.T[2] 

=#


function fst(alpha, liftfit, dcdalpha, alpha0)
    
    f =  (2*sqrt(abs(liftfit(alpha)/(dcdalpha*(alpha-alpha0)))) - 1)^2
    if f>= 1 || isnan(f)
        return 1.0
    elseif f<0
        return 0.0
    else
        return f
    end
    
end

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




function riso_state_rates(X, U, Udot, alpha, alphadot, c, dcldalpha, alpha0, afm, afp, liftfit, A1, A2, b1, b2, Tp, Tf) 

    ### Calculate constants
    Tu = c/(2*U) #Todo: This will go to NaN if U=0
    # println(typeof(X[1]))
    # @show U, Tu, c
    # @show X[1]
    # @show Tp
    
    ### Calculate the state rates
    dx1 = (b1*A1*alpha/Tu) - X[1]*(b1+ (c*Udot/(2*(U)^2)))/Tu 

    dx2 = (b2*A2*alpha/Tu) - X[2]*(b2+ (c*Udot/(2*(U)^2)))/Tu

    ae = alpha*(1-A1-A2) + X[1] + X[2] #Effective Angle of Attack
    # @show dcldalpha, ae, alpha0, Tu, alphadot, Tp, X[3]
    dx3 = (dcldalpha*(ae-alpha0) + pi*Tu*alphadot)/Tp - X[3]/Tp #TODO: Should this be able to go negative? 

    alphaf = (X[3]/dcldalpha)+alpha0 #Seperation Angle of Attack
    fp = fst(alphaf, liftfit, dcldalpha, alpha0) #Todo: Should I be using the fst function or the seperationpoint function? 
    dx4 = fp/Tf - X[4]/Tf 

    # @show dx1, dx2, dx3, dx4

    return SVector(dx1, dx2, dx3, dx4)
end




#Todo: I need to add the checks and balances from above. 
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





function riso_residual(dx, x, y, p, t, airfoil)  

    ### Extract inputs
    u, udot, v, vdot, theta, thetadot = y

    # println(y)

    ### Extract parameters
    c = p 

    ### Calculate Aero state rates and factors for BEM 
    d_X = riso_state_rates(x, u, udot, v, vdot, theta, thetadot, c, airfoil)

    ### Dynamic Stall Model State Rate Residuals
    r1 = dx[1] - d_X[1]
    r2 = dx[2] - d_X[2]
    r3 = dx[3] - d_X[3]
    r4 = dx[4] - d_X[4] 
    
    if isnan(r1)
        println("riso residual: ", [r1, r2, r3, r4])
        error("Nan in Riso")
    elseif isnan(r2)
        println("riso residual: ", [r1, r2, r3, r4])
        error("Nan in Riso")
    elseif isnan(r3)
        println("riso residual: ", [r1, r2, r3, r4])
        error("Nan in Riso")
    elseif isnan(r4)
        println("riso residual: ", [r1, r2, r3, r4])
        error("Nan in Riso")
    end
    # println("riso residual: ", [r1, r2, r3, r4])
    return SVector(r1, r2, r3, r4)
end


function get_riso_y(twist, env, frequency, amplitude, t)  

    u = env.Vinf(t)
    udot = env.Vinfdot(t)

    v = 0.0
    vdot = 0.0

    theta = twist + amplitude*cos(frequency*t)
    thetadot = -amplitude*frequency*sin(frequency*t)

    return [u, udot, v, vdot, theta, thetadot]
end

#Todo: Maybe make a riso ode function that takes in a p for a given time step. So with the intent of augmenting p every iteration. 


# Out-of-place ODE dispatch
function (model::Riso)(x, p, t)
    return riso_ode(model.detype, model, x, p, t)
end

# function riso_ode(detype::Iterative, x, p, t)
#     u, udot, v, vdot, theta, thetadot, c, dcldalpha, alpha0, afm, afp, liftfit, A1, A2, b1, b2, Tp, Tf = p

#     return riso_state_rates(x, u, udot, v, vdot, theta, thetadot, c, dcldalpha, alpha0, afm, afp, liftfit, A1, A2, b1, b2, Tp, Tf)
# end

function riso_ode(detype::Iterative, model, x, p, t)
    n = model.n

    #p = [u, udot, alpha, alphadot, c]
    u = view(p, 1:n)
    udot = view(p, n+1:2n)
    alpha = view(p, 2n+1:3n)
    alphadot = view(p, 3n+1:4n)
    c = view(p, 4n+1:5n)


    dx = Array{eltype(x), 1}(undef, 4*model.n)
    
    for i = 1:model.n
        idx = 4*(i-1)

        # @show u[i]

        xs = view(x, idx+1:idx+4)

        dx[idx+1:idx+4] = riso_state_rates(xs, u[i], udot[i], alpha[i], alphadot[i], c[i], model.airfoils[i].dcldalpha, model.airfoils[i].alpha0, model.airfoils[i].alphasep[1], model.airfoils[i].alphasep[2], model.airfoils[i].cl, model.airfoils[i].A[1], model.airfoils[i].A[2], model.airfoils[i].b[1], model.airfoils[i].b[2], model.airfoils[i].T[1], model.airfoils[i].T[2])
    end

    return dx
end


# Inplace ODE dispatch
function (model::Riso)(dx, x, p, t)
    return riso_ode!(model.detype::DEType, model, dx, x, p, t)
end


# Inplace, functional ode form of the Riso Model. 
function riso_ode!(detype::Functional, model, dx, x, p, t) 
    U, Udot, alpha, alphadot = p[1:4]
    cvec = view(p, 5:model.n+4)
    for i = 1:model.n
        airfoil = model.airfoils[i]
        idx = (i-1)*4
        xs = view(x, idx+1:idx+4)
        dxs = view(dx, idx+1:idx+4)
        riso_state_rates!(dxs, xs, U(t), Udot(t), alpha(t), alphadot(t), cvec[i], airfoil.A[1], airfoil.A[2], airfoil.b[1], airfoil.b[2], airfoil.dcldalpha, airfoil.alpha0, airfoil.T[1], airfoil.T[2], airfoil.cl, airfoil.alphasep[1], airfoil.alphasep[2])
    end
end





function static_riso_dae!(resids, dx, x, p, t) #Todo: I'm not sure this is correct. 
    U, Udot, alpha, alphadot, c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp = p
    @show typeof(x)

    riso_state_rates!(dx, x, U(t), Udot(t), alpha(t), alphadot(t), c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp)

    resids[1] = dx[1]
    resids[2] = dx[1]
    resids[3] = dx[3]
    resids[4] = dx[4]
end










# function Clfs(alpha, clfit, dcldalpha, alpha0, afm, afp)
#     fst = seperationpoint(alpha, afm, afp, clfit, dcldalpha, alpha0)
    
#     clst = clfit(alpha)
#     # @show fst
#     # if fst>=1 #Todo. This is supposed to happen automagically. 
#     #     return clst/2
#     # else
#     #     return (clst - dcldalpha*(alpha - alpha0)*fst)/(1-fst)
#     # end  
#     return (clst - dcldalpha*(alpha - alpha0)*fst)/(1-fst)
# end





"""
    Clfs(alpha, clfit, dcldalpha, alpha0)

The fully seperated coefficient of lift. The function has the seperation point function inserted, and simplified. 

## Inputs
- alpha - Current angle of attack (radians)
- clfit - A spline object to return the static coefficient of lift as a function of alpha
- dcldalpha - The slope of the attached region of the coefficient of lift vs alpha curve (1/radians). 
- alpha0 - The zero lift angle of attack (radians)
"""
function Clfs(alpha, clfit, dcldalpha, alpha0) #Todo. Determine which of these three functions to keep. I think I reworked the math so this one doesn't have issues. -> I think this is the one with the math worked differently. 
    cl = clfit(alpha)
    a = alpha - alpha0
    am = a*dcldalpha
    term = am*sqrt(abs(cl/(am))) - 3*cl

    # @show am
    return -am*term/(4*cl)
end

# function Clfs(alpha, liftfit, dcldalpha, alpha0)
#     f = fst(alpha, liftfit, dcldalpha, alpha0)
#     Cl = 0.0
#     if f>= 1.0 
#         Cl =  liftfit(alpha)/2
#     else
#         Cl = (liftfit(alpha) - (dcldalpha*(alpha-alpha0))*f)/(1-f)  
#     end

#     if isnan(Cl)
#         Cl = liftfit(alpha)/2
#     end

#     return Cl
# end







function riso_coefficients(x, U, alpha, alphadot, c, airfoil::Airfoil)
    return riso_coefficients(x, U, alpha, alphadot, c, airfoil.cl, airfoil.cd, airfoil.dcldalpha, airfoil.alpha0, airfoil.A[1], airfoil.A[2], airfoil.alphasep[1], airfoil.alphasep[2])
end


function riso_coefficients(x, U, alpha, alphadot, c, clfit, cdfit, dcldalpha, alpha0, A1, A2, afm, afp)

    Cd0 = cdfit(alpha0)

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
    # fterm = (sqrt(fae)-sqrt(fpp))/2 - (fae-fpp)/4
    Cd_dyn = cdst + (alpha - alpha_E)*Cl_dyn + (cdst - Cd0)*(t1 - t2)
    return Cl_dyn, Cd_dyn
end






export parsesolution

function parsestates(dsmodel::Riso, xds, p)
    n = dsmodel.n

    u = view(p, 1:n) #Todo. I should just change u, and v to U... since that's what the model intakes. Why change about so much? -> But it allows me to have inflow and heaving. -> Why not just create a function that handles that? 
    # udot = view(p, n+1:2n)
    alpha = view(p, 2n+1:3n)
    alphadot = view(p, 3n+1:4n)
    c = view(p, 4n+1:5n)

    cl = zeros(dsmodel.n)
    cd = zeros(dsmodel.n)

    for i = 1:dsmodel.n

        idx = 4*(i-1)
        x_section = xds[idx+1:idx+4]
        # y = [u[i], udot[i], v[i], vdot[i], theta[i], thetadot[i]] 
        # @show y
        # @show c
        cl[i], cd[i] = riso_coefficients(x_section, u[i], alpha[i], alphadot[i], c[i], dsmodel.airfoils[i])
        # cl[i], cd[i] = riso_coefs(x_section, y, c[i], dsmodel.airfoils[i])
    end
    
    return cl, cd
end



function parsesolution(dsmodel::Riso, sol, p)

    ### Unpack
    x = Array(sol)'
    t = sol.t
    U, Udot, alpha, alphadot = p[1:4]

    cvec = view(p, 5:dsmodel.n+4)

    nt = length(t)

    clvec = zeros(nt, dsmodel.n)
    cdvec = zeros(nt, dsmodel.n)

    ### run through the time steps and calculate the dynamic lift and drag (based on the states)
    for i = 1:nt
        ti = t[i]
        for j = 1:dsmodel.n
            clvec[i,j], cdvec[i,j] = riso_coefficients(x[i,:], U(ti), alpha(ti), alphadot(ti), cvec[j], dsmodel.airfoils[j])
        end
    end

    return clvec, cdvec, t
end







### Includes moment. 
# function parsesolution(sol, p, polar)
#     u = reduce(hcat, sol.u)'
#     n,m = size(u)
    
#     U, Udot, alpha, alphadot, c, A[1], A[2], b[1], b[2], dcldalpha, alpha0, Tp, Tf, liftfit, afm, afp = p

#     dragfit = Akima(polar[:,1], polar[:,3])
#     momentfit = Akima(polar[:,1], polar[:,4])

#     ### Find a^st, the distance between the center of pressure and the quarter chord. #Todo: Use this to get coefficient of moment in riso_coefficients. 
#     alphavec = polar[:,1]
#     nn = length(alphavec)
#     fvec = zeros(nn)
#     astvec = zeros(nn)
#     for i=1:nn 
#         fvec[i] = seperationpoint(alphavec[i], afm, afp, liftfit, dcldalpha, alpha0)
#         astvec[i] = (momentfit(alphavec[i])-momentfit(alpha0))/liftfit(alphavec[i])
#     end

    
#     mat = hcat(reverse(fvec), reverse(astvec))
#     mat = uniquemat!(mat)
    
#     affit = Akima(mat[:,1], mat[:,2])


    
#     Cl = zeros(n)
#     Cd = zeros(n)
#     Cm = zeros(n)
    
#     for i = 1:n
#         t = sol.t[i]

#         a34 = alpha(t)
#         ae = a34*(1-A[1]-A[2]) + u[i,1] + u[i,2]
#         Tu = c/(2*U(t))
        
#         clfs = Clfs(ae, liftfit, dcldalpha, alpha0) 
        
#         Cl[i] = dcldalpha*(ae-alpha0)*u[i,4] + clfs*(1-u[i,4]) + pi*Tu*alphadot(t) 
#         fae = seperationpoint(ae, afm, afp, liftfit, dcldalpha, alpha0)
#         fterm = (sqrt(fae)-sqrt(u[i,4]))/2 - (fae-u[i,4])/4
#         Cd[i] = dragfit(ae) + (alpha(t)-ae)*Cl[i] + (dragfit(ae)-dragfit(alpha0))*fterm
#         aterm = affit(u[i,4]) - affit(seperationpoint(ae, afm, afp, liftfit, dcldalpha, alpha0))
#         # println(aterm)
#         Cm[i] = momentfit(ae) + Cl[i]*(aterm) - pi*Tu*alphadot(t)/2
#     end
#     t = sol.t
#     return Cl, Cd, Cm, u, t
# end









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






