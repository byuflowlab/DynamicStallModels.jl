#=
Assuming piecewise constant terms
=#

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

function take_step!(xj, xjm, deltat, uj, ujm, udotj, udotjm, alphaj, alphajm, alphadotj, alphadotjm, c, dcldalpha, alpha0, afm, afp, A1, A2, b1, b2, Tp, Tf, clfit)
    Pi = b1*((uj + ujm)/c) + (udotj + udotjm)/(uj + ujm)
    Qi = b1*A1*(ujm*alphajm + uj*alphaj)/c 
    Ci = exp(-Pi*deltat)
    Ii = Qi*(1-Ci)/Pi
    xj[1] = Ci*xjm[1] + Ii

    Pi = b2*((uj + ujm)/c) + (udotj + udotjm)/(uj + ujm)
    Qi = b2*A2*(ujm*alphajm + uj*alphaj)/c
    Ci = exp(-Pi*deltat)
    Ii = Qi*(1-Ci)/Pi
    xj[2] = Ci*xjm[2] + Ii

    alphaEjm = alphajm*(1-A1-A2) + xjm[1] + xjm[2]
    alphaEj = alphaj*(1-A1-A2) +xj[1] + xj[2]
    t1 = dcldalpha*(alphaEjm + alphaEj - 2*alpha0)
    t2 = pi*c*(alphadotjm/ujm + alphadotj/uj)/2

    Pi = 1/Tp
    Qi = (t1 + t2)/(2*Tp)
    Ci = exp(-Pi*deltat)
    Ii = Qi*(1-Ci)/Pi
    xj[3] = Ci*xjm[3] + Ii

    alphafj = xj[3]/dcldalpha + alpha0 
    alphafjm = xjm[3]/dcldalpha + alpha0 
    fj = seperationpoint(alphafj, afm, afp, clfit, dcldalpha, alpha0)
    fjm = seperationpoint(alphafjm, afm, afp, clfit, dcldalpha, alpha0)

    Pi = 1/Tf
    Qi = (fj + fjm)/(2*Tf)
    Ci = exp(-Pi*deltat)
    Ii = Qi*(1-Ci)/Pi
    xj[4] = Ci*xjm[4] + Ii
    return xj
end

function indicialsolve(tvec, u, udot, alpha, alphadot, c, dcldalpha, alpha0, afm, afp, A, b, Tp, Tf, clfit, x0)
    nt = length(tvec)
    states = zeros(nt, 4)
    states[1,:] = x0

    for i = 2:nt
        dt = tvec[i]-tvec[i-1]

        states[i,:] = take_step!(states[i,:], states[i-1,:], dt, u(tvec[i]), u(tvec[i-1]), udot(tvec[i]), udot(tvec[i-1]), alpha(tvec[i]), alpha(tvec[i-1]), alphadot(tvec[i]), alphadot(tvec[i-1]), c, dcldalpha, alpha0, afm, afp, A[1], A[2], b[1], b[2], Tp, Tf, clfit)
    end

    return states
end