function prepenvironment(; c=0.1, M=0.3, a=343.3, amp=10.0, shift=10.0, k=0.1)
    v = M*a
    omega = k*2*v/c
    
    function U(t)
        return M*a
    end

    function Udot(t)
        return 0.0
    end

    function alpha(t)
        alf = shift + amp*sin(omega*t)
        return alf*(pi/180)
    end

    function alphadot(t)  
        alfd = amp*omega*cos(omega*t)
        return alfd*(pi/180)
    end

    function alphaddot(t)
        alfdd = -amp*(omega^2)*sin(omega*t)
        return alfdd*(pi/180)
    end
    
    return U, Udot, alpha, alphadot, alphaddot
end

function extractdata(sol)
    t = sol.t
    u = reduce(hcat, sol.u)'
    n,m = size(u)
    du = zeros(n,m)

    for j=1:m
        du[:,j] = gradient(t, u[:,j], t)
    end
    return t, u, du
end

"""
    myrk4(ydot, y0, rng, p;numsteps=100)

A simple implementation of the Runge-Kutta 4th order solver. 

### Inputs:
- ydot - Function representing state rate equations. The input to this function must be of form (dx, x, p, t) and must be in-place. 
- y0 - Vector of initial conditions
- rng - Range to integrate the functions across
- p - Parameters of the ODEs to be solved. Can be float, vector, object or whatever. I'm not calculating sensitivities, so it can also be time varying.
- numsteps - Integer - Number of time steps to take. It just does a linear distribution of time steps... nothing fancy. 

### Outputs:
- t - Vector of floats - time at which the ODEs were solved
- y - Array of floats - Column by column solution of each state variable across time (time is down, state variables left to right)
"""
function myrk4(ydot, y0, rng, p;numsteps=100, dy0 = [])
    n = length(y0)
    y = zeros(numsteps, n)
    dy = zeros(n)
    dY = zeros(numsteps, n)
    if length(dy0)>0
        dy[:] = dy0
    end
    # println(dy)
    

    # deltat = (rng[2]-rng[1])/numsteps
    # tn = rng[1]
    t = collect(range(rng[1], rng[2], length=numsteps))
    deltat = t[2]-t[1]

    y[1,:] = y0
    ydot(dy, y[1,:], p, t[1])
    dY[1,:] = deepcopy(dy)
    for i = 2:numsteps
        ydot(dy, y[i-1,:], p, t[i-1])
        k1 = deltat.*dy
        ydot(dy, y[i-1,:] .+ (k1./2), p, t[i-1]+(deltat/2))
        k2 = deltat.*dy
        ydot(dy, y[i-1, :] .+ (k2./2), p, t[i-1]+(deltat/2))
        k3 = deltat.*dy
        ydot(dy, y[i-1, :] .+ k3, p, t[i-1] + deltat)
        k4 = deltat.*dy

        y[i,:] = y[i-1, :] .+ (k1 .+ 2.0.*k2 .+ 2.0.*k3 .+ k4 )./6
        ydot(dy, y[i,:], p, t[i]) #I want to get the derivative at the current point. but I'm pretty sure that the BEM uses the derivative value.... so I don't know if the input value would be problematic... Also... I don't know if I should use y value at the current step... I think so. 
        # println(dy)
        dY[i,:] = deepcopy(dy)
    end
    return t, y, dY
end