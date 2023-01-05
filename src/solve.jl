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

        states[1,idxs], loads[1,:], p = initialize(dsmodel, Uvec, aoavec, tvec, airfoil, c, a)

        for i = 1:nt-1 
            ### Update environmental inputs
            t = tvec[i]
            dt = tvec[i+1]-tvec[i]
            # p[23] = U = Uvec[i+1]
            # p[24] = aoa = aoavec[i+1]
            update_environment!(dsmodel, p, Uvec[i+1], aoavec[i+1])

            if verbose
                @show t
            end
            
            if isa(dsmodel, BeddoesLeishman)
                if dsmodel.version==3
                    newstates = view(states, i+1, idxs)
                    dsmodel(states[i,idxs], newstates, p, t, dt)
                end
            else
                states[i+1, idxs] = dsmodel(states[i,idxs], p, t, dt)
            end

            loads[i+1,:] = getloads(dsmodel, states[i,idxs], p, airfoil)
        end
    end

    return states, loads
end