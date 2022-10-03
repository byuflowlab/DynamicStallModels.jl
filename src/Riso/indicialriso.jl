#=
Assuming piecewise constant terms
=#



function take_step!(xj, xjm, deltat, uj, ujm, udotj, udotjm, alphaj, alphajm, alphadotj, alphadotjm, c, dcldalpha, alpha0, afm, afp, A1, A2, b1, b2, Tp, Tf, clfit, airfoil) #TODO: Need to figure out life with the airfoil and now with repeated inputs. 
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
    # fj = seperationpoint(alphafj, afm, afp, clfit, dcldalpha, alpha0)
    # fjm = seperationpoint(alphafjm, afm, afp, clfit, dcldalpha, alpha0)
    fj = separationpoint(airfoil, alphafj)
    fjm = separationpoint(airfoil, alphafjm)
    

    Pi = 1/Tf
    Qi = (fj + fjm)/(2*Tf)
    Ci = exp(-Pi*deltat)
    Ii = Qi*(1-Ci)/Pi
    xj[4] = Ci*xjm[4] + Ii
    return xj
end

function indicialsolve(tvec, u, udot, alpha, alphadot, c, dcldalpha, alpha0, afm, afp, A, b, Tp, Tf, clfit, x0, airfoil)
    nt = length(tvec)
    states = zeros(nt, 4)
    states[1,:] = x0

    for i = 2:nt
        dt = tvec[i]-tvec[i-1]

        states[i,:] = take_step!(states[i,:], states[i-1,:], dt, u(tvec[i]), u(tvec[i-1]), udot(tvec[i]), udot(tvec[i-1]), alpha(tvec[i]), alpha(tvec[i-1]), alphadot(tvec[i]), alphadot(tvec[i-1]), c, dcldalpha, alpha0, afm, afp, A[1], A[2], b[1], b[2], Tp, Tf, clfit, airfoil)
    end

    return states
end

function parseindicialsolution(u, p, polar, airfoil)
    n,m = size(u)
    
    U, Udot, alpha, alphadot, c, A[1], A[2], b[1], b[2], dcldalpha, alpha0, Tp, Tf, liftfit, afm, afp = p

    dragfit = Akima(polar[:,1], polar[:,3])
    momentfit = Akima(polar[:,1], polar[:,4])

    ### Find a^st, the distance between the center of pressure and the quarter chord. #Todo: Use this to get coefficient of moment in riso_coefficients. 
    alphavec = polar[:,1]
    nn = length(alphavec)
    fvec = zeros(nn)
    astvec = zeros(nn)
    for i=1:nn 
        # fvec[i] = seperationpoint(alphavec[i], afm, afp, liftfit, dcldalpha, alpha0)
        fvec[i] = separationpoint(airfoil, alphavec[i])
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
        # fae = seperationpoint(ae, afm, afp, liftfit, dcldalpha, alpha0)
        fae = separationpoint(airfoil, ae)
        fterm = (sqrt(fae)-sqrt(u[i,4]))/2 - (fae-u[i,4])/4
        Cd[i] = dragfit(ae) + (alpha(t)-ae)*Cl[i] + (dragfit(ae)-dragfit(alpha0))*fterm
        # aterm = affit(u[i,4]) - affit(seperationpoint(ae, afm, afp, liftfit, dcldalpha, alpha0))
        aterm = affit(u[i,4]) - affit(separationpoint(airfoil, ae)) #TODO: Isn't the call to separationpoint() redundant because shouldn't that be the same thing as fae? 
        # println(aterm)
        Cm[i] = momentfit(ae) + Cl[i]*(aterm) - pi*Tu*alphadot(t)/2
    end
    return Cl, Cd, Cm
end