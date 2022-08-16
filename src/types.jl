export simpleairfoil,  complexairfoil, riso

struct Airfoil{TF, Tfit}
    polar::Array{TF, 2}
    cl::Tfit
    cd::Tfit
    cm::Tfit
    dcldalpha::TF
    alpha0::TF
    alphasep::Array{TF,1}
    A::Array{TF,1}
    b::Array{TF,1}
    T::Array{TF,1}
end


function simpleairfoil(polar)
    cl = Akima(polar[:,1], polar[:,2])
    cd = Akima(polar[:,1], polar[:,3])
    cm = Akima(polar[:,1], zeros(length(polar[:,1])))
    dcldalpha = 2*pi
    alpha0 = 0.0
    A = [0.3, 0.7] #Todo: Maybe do tuples... Since they aren't going to change... That might make it faster. 
    b = [0.14, 0.53]
    T = [1.7, 3.0]

    _, minclidx = findmin(polar[:,2]) #Todo: I have the find_seperation_alpha function. 
    _, maxclidx = findmax(polar[:,2])

    alphasep = [polar[minclidx, 1], polar[maxclidx,1]]
    return Airfoil(polar, cl, cd, cm, dcldalpha, alpha0, alphasep, A, b, T)
end

function nearestto(xvec, x) #TODO: Move to a utilities file. 
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

function complexairfoil(polar; A = [0.3, 0.7], b = [0.14, 0.53], T = [1.7, 3.0])
    #Todo: Need some sort of behavior when the provided polar is too small. 

    cl = Akima(polar[:,1], polar[:,2])
    cd = Akima(polar[:,1], polar[:,3])
    cm = Akima(polar[:,1], zeros(length(polar[:,1])))
    alpha0, _ = brent(cl, -0.25, 0.25)
    _, minclidx = findmin(polar[:,2])
    _, maxclidx = findmax(polar[:,2])

    alphasep = [polar[minclidx, 1], polar[maxclidx,1]]

    # @show minclidx, maxclidx
    middlepolar = polar[minclidx:maxclidx,:]
    # @show middlepolar
    _, cl0idx = nearestto(middlepolar[:,2], 0.0)
    alpha50 = middlepolar[end,1]*0.25

    # @show alpha50
    _, alf50idx = nearestto(middlepolar[:,1], alpha50)

    # @show cl0idx, alf50idx

    _, dcldalpha = linear_fit(middlepolar[cl0idx:alf50idx,1], middlepolar[cl0idx:alf50idx,2]) #TODO: Create my own linear fit function so I don't have to pull in a package. #Todo: This is returning a NaN
    if isnan(dcldalpha)
        dcldalpha=2*pi
        @warn("dcldalpha returned NaN")
    end
    return Airfoil(polar, cl, cd, cm, dcldalpha, alpha0, alphasep, A, b, T)
end





abstract type DEType end #Is the model designed to be solved in one go, or iteratively, updating p ever iteration. 

struct Functional <: DEType
end

struct Iterative <: DEType
end




abstract type DSModel end

struct Riso <: DSModel
    detype::DEType 
    n::Int #Number of airfoils simulated
    airfoils::Array{Airfoil,1}
end

function riso(airfoils) #::Array{Airfoil,1} #Todo: I'm not sure how to type this. 
    n = length(airfoils)
    return Riso(Iterative(), n, airfoils)
end

struct IndicialRiso <: DSModel
end

struct BeddoesLeishman <: DSModel
end

struct IndicialBeddoesLeishman <: DSModel
end