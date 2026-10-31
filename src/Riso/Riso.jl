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

export Riso
"""
    Riso(airfoils; detype=Continuous(), A, b, T)

Aggregate continuous Risø rate wrapper. A and b contain the two wake-response
amplitudes and decay coefficients; T contains pressure/separation times in
seconds. Constants are explicit, shared by these stations, and are not inferred
from unrelated Airfoil model fields. Each station has four states.

### Inputs
- detype - Continuous(); the aggregate wrapper is not a discrete state updater.
- airfoils - A vector of Airfoil structs, one corresponding to each section to be simulated. 
"""
struct Riso{TI,TA,TB,TT} <: DSModel
    detype::DEType 
    n::TI #Number of airfoils simulated
    airfoils::Array{Airfoil,1}
    A::TA
    b::TB
    T::TT
end

"""
    Riso(airfoils; detype=Continuous(), A, b, T)

A convenience. constructor for the Risø model. 

### Inputs
- airfoils - A vector of Airfoil structs, one corresponding to each section to be simulated. 
- A, b, T - Explicit two-entry vectors; T values are dimensional seconds.
"""
function Riso(airfoils; detype::DEType=Continuous(), A=nothing, b=nothing, T=nothing)
    detype isa Continuous || throw(ArgumentError("The aggregate Riso wrapper returns continuous rates; use Continuous()."))
    any(isnothing,(A,b,T)) && throw(ArgumentError("Supply explicit Riso A, b and T vectors; constants are not Airfoil fields. T contains Tp and Tf in seconds."))
    validate_riso_constants(A,b,T)
    n = length(airfoils)
    n>0 || throw(ArgumentError("Provide at least one airfoil."))
    return Riso(detype, n, Airfoil[airfoils...], copy(A), copy(b), copy(T))
end

function validate_riso_constants(A,b,T)
    all(length(v)==2 for v in (A,b,T)) || throw(DimensionMismatch("Riso A, b and T each need two entries."))
    all(v->all(isfinite,v),(A,b,T)) || throw(ArgumentError("Riso constants must be finite."))
    all(b .> 0) && all(T .> 0) || throw(ArgumentError("Riso decay constants and times must be positive."))
    return nothing
end

function riso_environment(y,t)
    env=evaluate_environment(y,t)
    all(isfinite,env) && env[1]>0 || throw(ArgumentError("Riso requires finite inputs and positive relative speed."))
    return env
end







function riso_state_rates(X, U, Udot, alpha, alphadot, c, dcldalpha, alpha0, afm, afp, liftfit, A1, A2, b1, b2, Tp, Tf, airfoil) #Passing the airfoil and the coefficients seems redundant. 

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
    # fp = fst_riso(alphaf, liftfit, dcldalpha, alpha0) #Todo: Should I be using the fst function or the separationpoint_riso function? 
    fp = separationpoint(airfoil, alphaf)
    dx4 = fp/Tf - X[4]/Tf 

    # @show dx1, dx2, dx3, dx4

    return SVector(dx1, dx2, dx3, dx4)
end




#Todo: I need to add the checks and balances from above. 
function riso_state_rates!(dx, x, U, Udot, alpha, alphadot, c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp, airfoil)

    # @show alpha
    Tu = c/(2*U)
    alpha_E = alpha*(1 - A1 - A2) + x[1] + x[2]
    alpha_f = x[3]/dcldalpha + alpha0
    # @show alpha_f
    # fst = separationpoint_riso(alpha_f, afm, afp, clfit, dcldalpha, alpha0)
    fst = separationpoint(airfoil, alpha_f)

    dx[1] = -x[1]*(b1 + c*Udot/(2*(U^2)))/Tu + b1*A1*alpha/Tu
    dx[2] = -x[2]*(b2 + c*Udot/(2*(U^2)))/Tu + b2*A2*alpha/Tu
    dx[3] = -x[3]/Tp + (dcldalpha*(alpha_E - alpha0) + pi*Tu*alphadot)/Tp
    dx[4] = -x[4]/Tf + fst/Tf
end





function riso_residual(dx, x, y, p, t, airfoil; A=nothing,b=nothing,T=nothing)
    length(dx)==length(x)==4 || throw(DimensionMismatch("Riso residual needs four states and four rates."))
    any(isnothing,(A,b,T)) && throw(ArgumentError("Supply explicit A, b and T for the Riso residual."))
    validate_riso_constants(A,b,T)
    if length(y)==6
        iszero(y[3]) && iszero(y[4]) || throw(ArgumentError("Convert heaving inputs to relative U/alpha before using the Riso residual."))
        y=[y[1],y[2],y[5],y[6]]
    end
    length(y)==4 || throw(DimensionMismatch("Riso environment is U, Udot, alpha, alphadot."))
    U,Udot,alpha,alphadot=riso_environment(y,t)
    c=p isa Number ? p : only(p)
    isfinite(c) && c>0 || throw(ArgumentError("Chord must be finite and positive."))
    rates=riso_state_rates(x,U,Udot,alpha,alphadot,c,airfoil.dcldalpha,airfoil.alpha0,airfoil.alphasep[1],airfoil.alphasep[2],airfoil.cl,A[1],A[2],b[1],b[2],T[1],T[2],airfoil)
    return SVector{4}(dx.-rates)
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
function (model::Riso)(x, p, t) #Todo: This is probably needing a revamp to x, y, p, t. -> Need to figure out how to do the dispatching. 
    dx=similar(x)
    model(dx,x,p,t)
    return dx
end


riso_ode(::Discrete,args...)=throw(ArgumentError("Use the Continuous() Riso rate wrapper."))

# Inplace ODE dispatch
function (model::Riso)(dx, x, p, t)
    return riso_ode!(model.detype::DEType, model, dx, x, p, t)
end


# Inplace, continuous ode form of the Riso Model. 
function riso_ode!(detype::Continuous, model, dx, x, p, t) 
    length(p)==4+model.n || throw(DimensionMismatch("Riso p contains four environment entries and one chord per airfoil."))
    length(dx)==length(x)==4model.n || throw(DimensionMismatch("Riso has four states per airfoil."))
    U, Udot, alpha, alphadot = riso_environment(p[1:4],t)
    cvec = view(p, 5:model.n+4)
    for i = 1:model.n
        airfoil = model.airfoils[i]
        idx = (i-1)*4
        xs = view(x, idx+1:idx+4)
        dxs = view(dx, idx+1:idx+4)
        isfinite(cvec[i]) && cvec[i]>0 || throw(ArgumentError("Chord must be finite and positive."))
        riso_state_rates!(dxs,xs,U,Udot,alpha,alphadot,cvec[i],model.A[1],model.A[2],model.b[1],model.b[2],airfoil.dcldalpha,airfoil.alpha0,model.T[1],model.T[2],airfoil.cl,airfoil.alphasep[1],airfoil.alphasep[2],airfoil)
    end
end





function static_riso_dae!(resids, dx, x, p, t, airfoil)
    U, Udot, alpha, alphadot, c, A1, A2, b1, b2, dcldalpha, alpha0, Tp, Tf, clfit, afm, afp = p
    U,Udot,alpha,alphadot=riso_environment([U,Udot,alpha,alphadot],t)
    rates=similar(x)
    riso_state_rates!(rates,x,U,Udot,alpha,alphadot,c,A1,A2,b1,b2,dcldalpha,alpha0,Tp,Tf,clfit,afm,afp,airfoil)
    resids .= dx .- rates
    return resids
end

static_riso_dae!(resids,dx,x,p,t)=throw(ArgumentError("Pass the Airfoil as the sixth argument to static_riso_dae!."))










# function Clfs(alpha, clfit, dcldalpha, alpha0, afm, afp)
#     fst = separationpoint_riso(alpha, afm, afp, clfit, dcldalpha, alpha0)
    
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
    # At simultaneous zero static/inviscid lift, the attached linear limit is zero.
    iszero(am) && iszero(cl) && return zero(am+cl)
    term = am*sqrt(abs(cl/(am))) - 3*cl

    # @show am
    return -am*term/(4*cl)
end

# function Clfs(alpha, liftfit, dcldalpha, alpha0)
#     f = fst_riso(alpha, liftfit, dcldalpha, alpha0)
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







function riso_coefficients(x, U, alpha, alphadot, c, airfoil::Airfoil; A=nothing)
    isnothing(A) && throw(ArgumentError("Supply the two Riso A coefficients explicitly."))
    return riso_coefficients(x,U,alpha,alphadot,c,airfoil.cl,airfoil.cd,airfoil.dcldalpha,airfoil.alpha0,A[1],A[2],airfoil.alphasep[1],airfoil.alphasep[2],airfoil)
end


function riso_coefficients(x, U, alpha, alphadot, c, clfit, cdfit, dcldalpha, alpha0, A1, A2, afm, afp, airfoil)

    Cd0 = cdfit(alpha0)

    Tu = c/(2*U)
    alpha_E = alpha*(1 - A1 - A2) + x[1] + x[2]

    # alpha_f = x[3]/dcldalpha + alpha0
    # fst = separationpoint_riso(alpha_f, afm, afp, clfit, dcldalpha, alpha0)

    ### Calculate Clfs #Todo: Should Clfs use fst calculated from alpha_f or alpha? I'm going to assume it is based off of alpha, but recognize that this could be wrong. 
    # cl_fs = Clfs(alpha_E, clfit, dcldalpha, alpha0, afm, afp)
    cl_fs = Clfs(alpha_E, clfit, dcldalpha, alpha0)
    # @show cl_fs

    Cl_dyn = dcldalpha*(alpha_E - alpha0)*x[4] + cl_fs*(1-x[4]) + pi*Tu*alphadot

    cdst = cdfit(alpha_E)
    # fst = separationpoint_riso(alpha_E, afm, afp, clfit, dcldalpha, alpha0)
    fst = separationpoint(airfoil, alpha_E)
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
        cl[i], cd[i] = riso_coefficients(x_section,u[i],alpha[i],alphadot[i],c[i],dsmodel.airfoils[i];A=dsmodel.A)
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
            u,_,a,adot=riso_environment([U,Udot,alpha,alphadot],ti)
            clvec[i,j], cdvec[i,j] = riso_coefficients(x[i,4(j-1)+1:4j],u,a,adot,cvec[j],dsmodel.airfoils[j];A=dsmodel.A)
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
#         fvec[i] = separationpoint_riso(alphavec[i], afm, afp, liftfit, dcldalpha, alpha0)
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
#         fae = separationpoint_riso(ae, afm, afp, liftfit, dcldalpha, alpha0)
#         fterm = (sqrt(fae)-sqrt(u[i,4]))/2 - (fae-u[i,4])/4
#         Cd[i] = dragfit(ae) + (alpha(t)-ae)*Cl[i] + (dragfit(ae)-dragfit(alpha0))*fterm
#         aterm = affit(u[i,4]) - affit(separationpoint_riso(ae, afm, afp, liftfit, dcldalpha, alpha0))
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
