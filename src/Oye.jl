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
"""
struct Oye{TI} <: DSModel
    detype::DEType 
    n::TI #Number of airfoils simulated
    airfoils::Array{Airfoil,1}
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
        
            c, dcndalpha, alpha0, A, U, aoa = ps #Inputs 
            
            xs = view(x, ns*(i-1)+1:ns*i) #Nodal states. 
            
            idx = ns*(i-1)+1:ns*i
            newstates[idx] = update_states(dsmodel, xs, U, aoa, dt, c, dcndalpha, alpha0, A, i)
        end
        return newstates
    end
end

function numberofstates(dsmodel::Oye)
    return 1
end

function numberofparams(dsmodel::Oye)
    return 6
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
        return getcoefficient_Oye_indicial(dsmodel, states, p, airfoil)
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
        envvars = [U, alpha]

        cn = airfoil.cn(alpha) #Static normal force
        cn_inv = airfoil.dcndalpha*(alpha-airfoil.alpha0) #Inviscid normal force 
        
        fst = (2*sqrt(abs(cn/cn_inv)) - 1)^2

        if fst>1
            fst = 1.0
        end
        
        states = [fst]

        A = airfoil.A[1]

        params = [c, airfoil.dcndalpha, airfoil.alpha0, A]

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

function update_states(dsmodel::Oye, oldstates, U, alpha, deltat, c, dcndalpha, alpha0, A, i)
    ### Unpack
    fold = oldstates[1]
    airfoil = dsmodel.airfoils[i]

    # @show alpha
    cn = airfoil.cn(alpha) #Static normal force
    cn_inv = dcndalpha*(alpha-alpha0) #Inviscid normal force 

    fst = (2*sqrt(abs(cn/cn_inv)) - 1)^2 #New separation point. #TODO: I wonder if I could hot swap this out for any separation point function of my choice. 
    tau = A*c/U #checked. 
    # @show tau

    f = fst + (fold - fst)*exp(-deltat/tau) #Delay on separation point #Todo: This isn't quite right. 

    if f>1
        return [1.0]
    else
        return [f]
    end
end


function getcoefficient_Oye_indicial(dsmodel::Oye, states, p, airfoil)
    ### Unpack
    f = states[1]
    _, dcndalpha, alpha0, _, _, alpha = p
    
    cn = airfoil.cn(alpha)
    cn_inv = dcndalpha*(alpha-alpha0)
    fst = (2*sqrt(abs(cn/cn_inv)) - 1)^2 
    cn_fs = (cn-(cn_inv*fst))/(1-fst) #Todo: I think this should actually use fst, not f. 

    Cn = f*cn_inv + (1-f)*cn_fs #Todo: Calculate the Cc, Cl, Cd, and potentially the Cm if possible. 

    return [Cn]
end

function updateenvironment_oye_indicial!(p, U, aoa)
    p[5] = U
    p[6] = aoa
end






