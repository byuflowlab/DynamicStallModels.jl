


function states!(dx, x, p, t)
    alpha, alphadot, liftfit, liftfitderiv, dcldalpha, A, omega, zeta = p

    dx[1] = -omega[1]*x[1] + omega[1]*dcldalpha*(alpha(t)) + (1-A[1])*dcldalpha*alphadot(t)
    dx[2] = x[3]

    delCl = dcldalpha*alpha(t) - liftfit(alpha(t))
    delCldot = dcldalpha*alphadot(t) - liftfitderiv(alpha(t))*alphadot(t)

    dx[3] = -2*zeta*omega[2]*x[3] - (omega[2]^2)*(1 + (zeta^2))*x[2] - (omega[2]^2)*(1 + (zeta^2))*(delCl + A[2]*delCldot)
end

function parsesolution(sol, p)
    u = reduce(hcat, sol.u)'
    n,m = size(u)

    alpha, alphadot, liftfit, liftfitderiv, dcldalpha, A, omega, zeta = p

    t = sol.t

    Cl = zeros(n)

    for i=1:n
        Cl[i] = u[i,1] + u[i,2]
    end
    return Cl, t, u
end