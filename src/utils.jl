function numberofstates_total(dsmodel::DSModel)
    return dsmodel.n*numberofstates(dsmodel)
end

function numberofparams_total(dsmodel::DSModel)
    return dsmodel.n*numberofparams(dsmodel)
end

function numberofloads_total(dsmodel::DSModel)
    return dsmodel.n*numberofloads(dsmodel)
end

function nearestto(xvec, x) 
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end



function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end



function uniquemat!(mat;column=1)
    n, m = size(mat)
    if column>m
        error("uniquemat!: You must choose a column within the matrix.")
    end
    listofindices = Int[]
    for i=1:length(mat[:,1])
        value = mat[i,column]
        index = findfirst(x->x==value,mat[:,column])
        push!(listofindices, index)
    end
    unique!(listofindices)
    return mat[listofindices,:]
end

function Heavi(t)
    if t>=0.0
        return 1.0
    else
        return 0.0
    end
end

function diracdelta(n;tol=0.001)  
    if isapprox(n, 0.0, atol=tol)
        return typeof(n)(1)
    else
        return typeof(n)(0)
    end
end

function dirac(x; w=0.1) 
    #Quadratic Dirac Delta (It's continuous baby, so it should be caught by DifferentialEquations.)
    L = w/2
    return L/(pi*(x^2 + L^2))
end


export Linear #Todo: Add to FLOWMath. 

struct Linear
    x
    y

    function Linear(x,y)
        for i in eachindex(x)[2:end-1]
            if !(x[i-1]<x[i]<x[i+1])
                error("Linear(): x vector must be in non-repeating ascending order.")
            end
        end
        
        new(x,y)
    end
end



function (interp::Linear)(x)
    if x<interp.x[1]
        @warn("Linear(): Outside of linear interpolation domain.")
        return interp.y[1]
    elseif x>interp.x[end]
        @warn("Linear(): Outside of linear interpolation domain.")
        return interp.y[end]
    end

    idx = findfirst(i -> x<i, interp.x)


    x0 = interp.x[idx-1]
    x1 = interp.x[idx]
    y0 = interp.y[idx-1]
    y1 = interp.y[idx]

    top = y0*(x1-x) + y1*(x-x0)
    bot = x1-x0

    return top/bot
end


