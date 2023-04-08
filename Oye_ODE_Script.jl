using DifferentialEquations, Plots


"""
M = Mach Number
a = Speed of Sound
A = vector for passed in constants
c = chord length
"""


function Oye_State_Rate!(dx, x, A, c, airfoil, Uvec, alphavec, t)
    T_f = (Uvec(t))/(A[1]*c)

    fst = separationpoint(airfoil, alphavec(t))

    dx[1] = -T_f*(x[1]-fst)
end

function Oye_ODE!(model, dx, x, p, t)
    A, c, Uvec, alphavec = p
    for i in 1:model.n
        airfoil = model.airfoils[i]
        idx = (i-1)*1
        xs = view(x , idx+1:idx+1)
        dxs = view(dx, idx+1:idx+1)
        Oye_State_Rate!(dxs, xs, A, c, airfoil, Uvec, alphavec, t)
    end
end

function (model::Oye)(dx, x , p, t)
    Oye_ODE!(model, dx, x, p, t)
end

function parsesolution_Oye(sol, p, af)
    _, _, _, alphavec = p
    f = Array(sol)
    tvec = sol.t
    airfoil = af

    alpha0 = airfoil.alpha0


    
    alpha = []
    Lift = []

    for i in 1:length(tvec)

        push!(alpha, alphavec(tvec[i]))
        
        C_inv = airfoil.dcldalpha*(alphavec(tvec[i]) - alpha0)

        cl = airfoil.cl(alphavec(tvec[i]))

        cl_sep = airfoil.cl(airfoil.alphasep[2])

        C_fs = cl_fullysep_faber(cl, cl_sep, airfoil.dcldalpha, alphavec(tvec[i]), alpha0, airfoil.alphasep[2])

        C_L_dyn = f[i]*C_inv+(1-f[i])*C_fs

        push!(Lift, C_L_dyn)
    end

    return Lift, alpha

end





#practice ODE
"""
function U(t)
    return 60.0
end

function Alpha(t)
    return 10*sin(6*t) + 8
end

p = [4, 0.55, 0.5]

function Oye_Practice(dx, x , p , t)
    T_f = U(t)/(p[1]*p[2])

    fst = p[3]

    dx[1] = -T_f*(x[1]-fst)
end


tspan = (0.0, 10.0)

x0 = [0.7]

prob = ODEProblem(Oye_Practice, x0, tspan, p)

sol = solve(prob)

plot(sol)
"""


