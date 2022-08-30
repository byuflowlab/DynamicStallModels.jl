

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
