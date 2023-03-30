#=
The Oye model, as given by ???. 

Indicial States
1 - f
=#
export Oye

"""
    Oye(detype::DEType, n::Int, airfoils::Array{Airfoil, 1})

The Øye model struct. It stores airfoil data for every section to be simulated. It can be used as a method to return updated states or state rates depending on it's DEType. 

### Inputs
- detype - The type of model it is, Functional(), Iterative(), or Indicial().
- n - The number of sections to be simulated. 
- airfoils - A vector of Airfoil structs, one corresponding to each section to be simulated.
- cflag::Int - A flag to apply the separation delay to the coefficient of 1) lift, 2) normal force. 
- version::Int - A flag to say whether to use 1) Hansen 2008, or 2) Faber 2018's implementation of the model.
"""
struct Oye{TI} <: DSModel
    detype::DEType 
    n::TI #Number of airfoils simulated
    airfoils::Array{Airfoil,1}
    cflag::TI 
    version::TI
end

function Oye(detype::DEType, airfoils; cflag::Int=1, version::Int=1)
    return Oye(detype, length(airfoils), airfoils, cflag, version)
end

function (dsmodel::Oye)(x, p, t, dt)
    if isa(dsmodel.detype, Functional)
        error("The state space Oye model is not setup yet.")
    elseif isa(dsmodel.detype, Iterative)
        error("The iterative Oye model is not set up yet.")
    else #The model is indicial
        nst = numberofstates_total(dsmodel)
        ns = numberofstates(dsmodel)
        np = numberofparams(dsmodel)
        newstates = Array{eltype(p), 1}(undef, nst)
        for i = 1:dsmodel.n
            ps = view(p, np*(i-1)+1:np*i)
        
            c, dcndalpha, alpha0, alphasep, A, U, aoa = ps #Inputs 
            
            xs = view(x, ns*(i-1)+1:ns*i) #Nodal states. 
            
            idx = ns*(i-1)+1:ns*i
            newstates[idx] = update_states(dsmodel, xs, U, aoa, dt, c, dcndalpha, alpha0, alphasep, A, i)
        end
        return newstates
    end
end

function numberofstates(dsmodel::Oye)
    return 1
end

function numberofparams(dsmodel::Oye)
    return 7
end

function numberofloads(dsmodel::Oye)
    return 1
end

function getloads(dsmodel::Oye, states, p, airfoil)
    if isa(dsmodel.detype, Functional)
        error("Oye functional implementation not yet prepared.")
    elseif isa(dsmodel.detype, Iterative)
        error("Oye iterative implementation not prepared for use yet.")
    else #Indicial
        if dsmodel.version == 1
            # println("Got here")
            return getcoefficient_indicial_hansen(dsmodel, states, p, airfoil)
        elseif dsmodel.version == 2
            return getcoefficient_indicial_faber(dsmodel, states, p, airfoil)
        else
            ver = dsmodel.version
            return getcoefficient_indicial_hansen(dsmodel, states, p, airfoil)
        end
    end
end

function initialize(dsmodel::Oye, Uvec, aoavec, tvec, airfoil::Airfoil, c, a)
    if isa(dsmodel.detype, Functional)
        @warn("Oye Functional implementation isn't prepared yet. - initialize()")
    elseif isa(dsmodel.detype, Iterative)
        @warn("Oye Iterative implementation isn't prepared yet. - initialize()")
    else #Model is indicial
        U = Uvec[1]
        alpha = aoavec[1]
        # @show alpha
        alphasep = airfoil.alphasep[2]
        alpha0 = airfoil.alpha0
        envvars = [U, alpha]

        if dsmodel.cflag == 2
            cn = airfoil.cn(alpha) #Static normal force
            cn_sep = airfoil.cn(alphasep)
            dcndalpha = airfoil.dcndalpha
        else
            cn = airfoil.cl(alpha)
            cn_sep = airfoil.cl(alphasep)
            dcndalpha = airfoil.dcldalpha
        end

        cn_inv = dcndalpha*(alpha-airfoil.alpha0) #Inviscid normal force 
        
        if dsmodel.version == 1
            fst = (2*sqrt(abs(cn/cn_inv)) - 1)^2
        elseif dsmodel.version == 2 #Todo: This will get eliminated here with implementation of multiple dispatch of separation point function. 
            cn_fs = cl_fullysep_faber(cn, cn_sep, dcndalpha, alpha, alpha0, alphasep)
            fst = (cn - cn_fs)/(cn_inv - cn_fs)
        else
            @warn("$ver is not an available version, defaulting to Hansen 2008 (Option 1).")
            fst = (2*sqrt(abs(cn/cn_inv)) - 1)^2
        end

        if fst>1
            fst = 1.0
        end
        
        states = [fst]

        

        A = airfoil.A[1]

        params = [c, dcndalpha, alpha0, alphasep, A]

        p = vcat(params, envvars)
        
        loads = getloads(dsmodel, states, p, airfoil) #Todo: Calculate Cld
         
        
        return states, loads, p
    end
end

function update_environment!(dsmodel::Oye, p, U, aoa)
    if isa(dsmodel.detype, Functional)
        @warn("Oye functional implementation isn't prepared yet. - initialize()")
    elseif isa(dsmodel.detype, Iterative)
        @warn("Oye iterative implementation isn't prepared yet. - initialize()")
    else #Model is indicial
        updateenvironment_oye_indicial!(p, U, aoa)
    end
end








######## Indicial code

function update_states(dsmodel::Oye, oldstates, U, alpha, deltat, c, dcndalpha, alpha0, alphasep, A, i)
    if dsmodel.version == 1
        return update_states_oye_hansen(dsmodel, oldstates, U, alpha, deltat, c, dcndalpha, alpha0, A, i)
    elseif dsmodel.version == 2 
        return update_states_oye_faber(dsmodel, oldstates, U, alpha, deltat, c, dcndalpha, alpha0, alphasep, A, i)
    else
        ver = dsmodel.version
        @warn("$ver not an option, defaulting to Hansen 2008 (option 1)")
    end
end

#=
Use Hansen's implementation of the Øye model (from his 2008 textbook)
=#
function update_states_oye_hansen(dsmodel::Oye, oldstates, U, alpha, deltat, c, dcndalpha, alpha0, A, i)
    ### Unpack
    fold = oldstates[1]
    airfoil = dsmodel.airfoils[i]

    # @show alpha
    if dsmodel.cflag == 2
        cn = airfoil.cn(alpha) #Static normal force
    else
        cn = airfoil.cl(alpha)
    end
    cn_inv = dcndalpha*(alpha-alpha0) #Inviscid normal force

    fst = (2*sqrt(abs(cn/cn_inv)) - 1)^2 #New separation point. #TODO: I wonder if I could hot swap this out for any separation point function of my choice. 
    tau = A*c/U #checked. 
    # @show tau

    f = fst + (fold - fst)*exp(-deltat/tau) #Delay on separation point #Todo: This isn't quite right. -> Is it wrong? It seems okay. 

    if f>1
        return [1.0]
    else
        return [f]
    end
end

#=
Use Faber's implementation of the Øye model (from his 2018 thesis), which is incidently also the Larsen 2007 implementation. 
=#
function update_states_oye_faber(dsmodel::Oye, oldstates, U, alpha, deltat, c, dcndalpha, alpha0, alphasep, A, i) #Todo: 
    ### Unpack
    fold = oldstates[1]
    airfoil = dsmodel.airfoils[i]

    # @show alpha
    if dsmodel.cflag == 2
        cn = airfoil.cn(alpha) #Static normal force
        cn_sep = airfoil.cn(alphasep)
    else
        cn = airfoil.cl(alpha)
        cn_sep = airfoil.cl(alphasep)
    end
    cn_inv = dcndalpha*(alpha-alpha0) #Inviscid normal force

    cn_fs = cl_fullysep_faber(cn, cn_sep, dcndalpha, alpha, alpha0, alphasep)
    fst = (cn - cn_fs)/(cn_inv - cn_fs)
    if alpha>alphasep
        fst = 0.0
    end
    # fst = separationpoint(airfoil, alpha)

    tau = A*c/U 

    f = fst + (fold - fst)*exp(-deltat/tau) #Delay on separation point

    if f>1
        return [1.0]
    elseif f<0
        return [0.0]
    else
        return [f]
    end
end

#=
Faber's fully separated lift coefficient from his 2018 thesis. 
=#
function cl_fullysep_faber(cl, cl_sep, dcldalpha, alpha, alpha0, alpha_sep)
    if alpha0 < alpha < alpha_sep
        #Hermite Interpolation as in Faber 2018
        #TODO: Needs to be extended to work for angles less than alpha0, it should have a reflection to what is done here. 
        t0 = (alpha - alpha0)/(alpha_sep - alpha0) #Faber 2018 EQ A.1a
        #As alpha -> alpha_sep this will diverge. So... what keeps it from diverging? 
        t1 = (alpha - alpha_sep)/(alpha_sep - alpha0) #Faber 2018 EQ A.1b
        term1 = (alpha_sep - alpha0)*dcldalpha*(1 + t0*(7*t1/6 - 1))/2
        term2 = cl_sep*t0*(1-( 2*t1))
        # return t0*(term1+term2) #Faber 2018 EQ A.4
        clfs = t0*(term1+term2) #Faber 2018 EQ A.4

        # if clfs >10
        #     @show t0, t1, term1, term2
        # end

        if t0>10
            @show alpha, alpha0, alpha_sep
        end

        return clfs
    else
        return cl
    end
end

function getcoefficient_indicial_hansen(dsmodel::Oye, states, p, airfoil)
    ### Unpack
    f = states[1]
    # c, dcndalpha, alpha0, alphasep, A, U, aoa = ps
    _, dcndalpha, alpha0, _, _, _, alpha = p
    
    if dsmodel.cflag == 2
        cn = airfoil.cn(alpha)
    else
        cn = airfoil.cl(alpha)
    end

    cn_inv = dcndalpha*(alpha-alpha0) #Hansen 2004 EQ 19.
    fst = (2*sqrt(abs(cn/cn_inv)) - 1)^2 #Hansen 2004 EQ 15.
    cn_fs = (cn-(cn_inv*fst))/(1-fst) #Hansen 2004 EQ 18.  

    Cn = f*cn_inv + (1-f)*cn_fs #Hansen 2004 EQ 17.
    #Todo: Calculate the Cc, Cl, Cd, and potentially the Cm if possible. 
    # @show dcndalpha, alpha, alpha0, f

    return [Cn]
end

function getcoefficient_indicial_faber(dsmodel::Oye, states, p, airfoil)
    ### Unpack
    f = states[1]
    # c, dcndalpha, alpha0, alphasep, A, U, aoa = ps
    _, dcndalpha, alpha0, alphasep, _, _, alpha = p
    
    if dsmodel.cflag == 2
        cn = airfoil.cn(alpha)
        cn_sep = airfoil.cn(alphasep)
    else
        cn = airfoil.cl(alpha)
        cn_sep = airfoil.cl(alphasep)
    end

    cn_inv = dcndalpha*(alpha-alpha0) #Hansen 2004 EQ 19.



    cn_fs = cl_fullysep_faber(cn, cn_sep, dcndalpha, alpha, alpha0, alphasep) #Faber 2018 EQ A.4

    Cn = f*cn_inv + (1-f)*cn_fs #Hansen 2004 EQ 17.
    #Todo: Calculate the Cc, Cl, Cd, and potentially the Cm if possible. 
    # if Cn>10
    #     # @show cn_inv, cn_fs, f
    #     @show cn, cn_sep, alpha
    # end

    return [Cn]
end

function updateenvironment_oye_indicial!(p, U, aoa)
    p[6] = U
    p[7] = aoa
end






