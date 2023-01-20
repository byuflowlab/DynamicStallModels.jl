#=
A place to hold all the solvers that aren't provided by DifferentialEquations.

Adam Cardoza 8/30/22
=#

export solve_indicial

function solve_indicial(dsmodel::DSModel, cvec, tvec, Uvec, aoavec; a = 343.0, verbose::Bool=false)

    nt = length(tvec)

    ### Initialize the states

    ns = numberofstates_total(dsmodel)
    states = Array{eltype(cvec), 2}(undef, nt, ns)
    nl = numberofloads_total(dsmodel)
    loads = Array{eltype(cvec), 2}(undef, nt, nl)

    # @show size(states), size(loads)
    
    for j = 1:dsmodel.n #Loop over each airfoil, n is the number of airfoils
        c = cvec[j]
        airfoil = dsmodel.airfoils[j]

        idxs = ns*(j-1)+1:j*ns 
        # @show idxs

        avec, bvec, p = initialize(dsmodel, Uvec, aoavec, tvec, airfoil, c, a) #! as of Jan 20, 2023 this needs to be updated

        # @show avec, bvec, p

        states[1,idxs] = avec #! as of Jan 20, 2023 this needs to be updated
        loads[1,:] = bvec

        # @show size(states), size(loads)

        for i = 1:nt-1 
            ### Update environmental inputs
            t = tvec[i]
            dt = tvec[i+1]-tvec[i]
            # p[23] = U = Uvec[i+1]
            # p[24] = aoa = aoavec[i+1]
            update_environment!(dsmodel, p, Uvec[i+1], aoavec[i+1])

            if verbose #telling it to be verbose will show what time step it is on
                @show t
            end
            
            if isa(dsmodel, BeddoesLeishman)
                if dsmodel.version==3
                    newstates = view(states, i+1, idxs)
                    dsmodel(states[i,idxs], newstates, p, t, dt) #in place function, this doesn't have to allocate memory every time step
                end
            else
                states[i+1, idxs] = dsmodel(states[i,idxs], p, t, dt) #solves for the new states, allocates memory every time step
            end

            loads[i+1,:] = getloads(dsmodel, states[i,idxs], p, airfoil)
        end
    end

    return states, loads
end