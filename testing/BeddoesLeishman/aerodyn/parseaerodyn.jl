
#Todo: I might just want to have one of these functions. 
function parseaerodyn(entries, items)
    ni = length(items)
    ne = Int(length(entries)/ni)
    mat = zeros(ne, ni)
    
    for i = 1:ne
        idx = (i-1)*ni+1:(i*ni)
        mat[i,:] = entries[idx]
    end
    return mat
end

