#=
A place to hold all the solvers that aren't provided by DifferentialEquations.

Adam Cardoza 8/30/22
=#

export solve_indicial

function solve_indicial(dsmodel::DSModel, cvec, tvec, Uvec, aoavec; a = 343.0, verbose::Bool=false)

    nt = length(tvec)

    ### Initialize the states
    ns = numberofstates(dsmodel)
    states = Array{eltype(cvec), 2}(undef, nt, ns)
    nl = numberofloads(dsmodel)
    loads = Array{eltype(cvec), 2}(undef, nt, nl)
    
    for j = 1:dsmodel.n
        c = cvec[j]
        airfoil = dsmodel.airfoils[j]

        idxs = ns*(j-1)+1:j*ns 

        states[1,idxs], loads[1,:], p = initialize(dsmodel, aoavec, tvec, airfoil, c, a)

        for i = 1:nt-1 
            
            if verbose
                t = tvec[i]
                @show t
            end
            dt = tvec[i+1]-tvec[i]
            U = Uvec[i+1]
            aoa = aoavec[i+1]

            ### Update environmental inputs
            y = [U, aoa, dt]

            states[i+1,idxs] = dsmodel(states[i,idxs], p, y) #update_states_ADO(dsmodel, states[i,idxs], flags, c, a, U, dt, aoa, dcndalpha, alpha0, A1, A2, b1, b2, Tf0, Tv0, Tp, Tvl, Cn1, alpha1, alpha2, S1, S2, S3, S4)

            loads[i+1,:] = getloads(dsmodel, states[i,idxs], p, y, airfoil)
        end
    end

    return states, loads
end