export Oye, parsesolution

"""
    Oye(detype::DEType, cflag::Int, version::Int, A::Real)

One-state Øye lift model: `Continuous()` rates or `Discrete()` exponential
steps. Use `cflag=1`; normal-force mode (`cflag=2`) is not implemented.
`version=1` selects Hansen's fully separated lift and `version=2` selects
Faber/Larsen's positive-incidence Hermite construction. The separation time is
`tau=A*c/U`, with positive A and chord, and nonnegative relative speed.

This lift-only model returns `(Cl, 0, 0)`, not predictions of drag or moment.
Angles are in radians. See the Øye examples for matching separation functions.
"""
struct Oye{TI, TF} <: DSModel
    detype::DEType
    cflag::TI
    version::TI
    A::TF
    function Oye{TI, TF}(detype::DEType, cflag::TI, version::TI, A::TF) where {TI, TF}
        cflag == 1 || throw(ArgumentError("Øye supports cflag=1 (lift); normal-force mode is not implemented."))
        version in (1, 2) || throw(ArgumentError("Øye version must be 1 (Hansen) or 2 (Faber/Larsen)."))
        isfinite(A) && A > 0 || throw(ArgumentError("Øye A must be finite and positive."))
        new{TI, TF}(detype, cflag, version, A)
    end
    Oye(detype::DEType, cflag::TI, version::TI, A::TF) where {TI, TF} = Oye{TI, TF}(detype, cflag, version, A)
end

numberofstates(::Oye) = 1
numberofparams(::Oye) = 2 # chord, xcp (unused by this lift-only model)
numberofloads(::Oye) = 3
get_cn(airfoil, alpha) = airfoil.model.cflag == 2 ? airfoil.cn(alpha) : airfoil.cl(alpha)
get_dcndalpha(airfoil) = airfoil.model.cflag == 2 ? airfoil.dcndalpha : airfoil.dcldalpha

function oye_attachment(airfoil, alpha)
    f = separationpoint(airfoil, alpha)
    isfinite(f) || throw(ArgumentError("Øye static attachment must be finite."))
    return clamp(f, zero(f), one(f))
end

function oye_rate(model::Oye, U, c)
    isfinite(U) && U >= 0 || throw(ArgumentError("Øye relative speed must be finite and nonnegative."))
    isfinite(c) && c > 0 || throw(ArgumentError("Øye chord must be finite and positive."))
    return U/(model.A*c) # 1/tau; the state freezes at U=0.
end

function initialize(model::Oye, airfoil::Airfoil, tvec, y, p)
    oye_rate(model, y[1], p[1])
    states = [oye_attachment(airfoil, y[3])]
    return states, collect(get_loads(model, airfoil, states, y, p))
end

function update_states(model::Oye, airfoil::Airfoil, oldstate, y, p, dt)
    newstate = similar(oldstate)
    return update_states!(model, airfoil, oldstate, newstate, y, p, dt)
end

function update_states!(model::Oye, airfoil::Airfoil, oldstate, newstate, y, p, dt)
    isfinite(dt) && dt >= 0 || throw(ArgumentError("Øye time step must be finite and nonnegative."))
    rate = oye_rate(model, y[1], p[1])
    fst = oye_attachment(airfoil, y[3])
    newstate[1] = fst + (oldstate[1]-fst)*exp(-dt*rate)
    return newstate
end

function get_loads(model::Oye, airfoil::Airfoil, states, y, p)
    alpha = y[3]
    cinv = get_dcndalpha(airfoil)*(alpha-airfoil.alpha0)
    cfs = model.version == 2 ? cl_fullysep_faber(airfoil, alpha) : cl_fullysep_hansen(airfoil, alpha)
    cl = states[1]*cinv + (1-states[1])*cfs
    return cl, zero(cl), zero(cl)
end

function get_loads!(model::Oye, airfoil::Airfoil, states, loads, y, p)
    loads .= get_loads(model, airfoil, states, y, p)
    return loads
end

function cl_fullysep_hansen(airfoil, alpha)
    cstatic = get_cn(airfoil, alpha)
    cinv = get_dcndalpha(airfoil)*(alpha-airfoil.alpha0)
    fst = oye_attachment(airfoil, alpha)
    # Removable 0/0: the Hansen relation's f->1 limit is Cinv/2.
    # For a custom f=1, Cfs is unidentifiable; use this same convention.
    # The fully attached load itself is independent of Cfs.
    isone(fst) && return cinv/2
    return (cstatic-cinv*fst)/(1-fst)
end

"""Faber/Larsen Hermite fully separated lift on alpha0 < alpha < alpha_sep.
Outside that positive-incidence interval, use the static polar. No mirrored
negative-stall calibration is introduced by this implementation.
"""
function cl_fullysep_faber(airfoil, alpha)
    alpha0 = airfoil.alpha0
    alpha_sep = airfoil.alphasep[2]
    if alpha0 < alpha < alpha_sep
        t0 = (alpha-alpha0)/(alpha_sep-alpha0)
        t1 = (alpha-alpha_sep)/(alpha_sep-alpha0)
        term1 = (alpha_sep-alpha0)*get_dcndalpha(airfoil)*(1+t0*(7*t1/6-1))/2
        term2 = get_cn(airfoil, alpha_sep)*t0*(1-2*t1)
        return t0*(term1+term2)
    end
    return get_cn(airfoil, alpha)
end

"""
    state_rates!(model::Oye, airfoil, dx, x, y, p, t)

Evaluate df/dt=(f_static(alpha(t))-f)/tau, tau=A*c/U(t).
The four entries of y are U, Udot, alpha, alphadot (numbers or functions).
p contains chord and xcp; derivatives and xcp are unused by Øye.
"""
function state_rates!(model::Oye, airfoil::Airfoil, dx, x, y, p, t)
    U, _, alpha, _ = evaluate_environment(y, t)
    dx[1] = oye_rate(model, U, p[1])*(oye_attachment(airfoil, alpha)-x[1])
    return dx
end

"""
    parsesolution(model::Oye, airfoil, sol, y, p)

Return two rows: angle in radians and dynamic Cl at sol.t. Uses the same
coefficient routine as the discrete solver. p contains chord and xcp.
"""
function parsesolution(model::Oye, airfoil::Airfoil, sol, y, p)
    states = Array(sol)
    result = zeros(eltype(states), 2, length(sol.t))
    for i in eachindex(sol.t)
        env = collect(evaluate_environment(y, sol.t[i]))
        result[1,i] = env[3]
        result[2,i] = get_loads(model, airfoil, view(states,:,i), env, p)[1]
    end
    return result
end
