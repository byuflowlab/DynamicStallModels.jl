

"""
pareseaerodyn(entries, items, numnodes)

For parsing intermediate results from AeroDyn, with multiple stations. 
"""
function parseaerodyn(entries, items, numnodes)

    ni = length(items)
    nt = Int(length(entries)/(ni*numnodes))

    if numnodes==1
        mat = zeros(nt, ni)
    
        for i = 1:nt
            idx = (i-1)*ni+1:(i*ni)
            mat[i,:] = entries[idx]
        end
        return mat
    else
        mat = zeros(nt, ni, numnodes)
        
        idx = 1
        for i = 1:nt
            for j=1:numnodes
                for k = 1:ni
                    mat[i,k,j] = entries[idx]
                    idx += 1
                end
            end
        end
        return mat
    end

    error("parseaerodyn(): Didn't parse the text file. ")
end

