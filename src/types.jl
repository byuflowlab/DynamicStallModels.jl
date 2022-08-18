export simpleairfoil, airfoil, riso

"""
    Airfoil(polar, cl, cd, cm, dcldalpha, alpha0, alphasep, A, b, T)

### Inputs
- polar::Array{TF, 2} - A matrix of floats describing the airfoil polar. This includes the angle of attack (radians), and coefficients of lift, drag, and moment. 
- cl::Tfit - A fit of the coefficient of lift as a function of the angle of attack (radians). 
- cd::Tfit - A fit of the coefficient of drag as a function of the angle of attack (radians). 
- cm::Tfit - A fit of the coefficient of moment as a function of the angle of attack (radians). 
- dcldalpha::TF - The lift curve slope in the linear region (1/radians). Typically near 2 pi. 
- alpha0::TF - The zero lift angle of attack (radians). 
- alphasep::Array{TF, 1} - A vector of floats holding the angles of attack at which flow fully seperates from the airfoil. In order from least to greatest. 
- A::Array{TF, 1} - A vector of floats holding the A dynamic constants for the airfoil. 
- b::Array{TF, 1} - A vector of floats holding the b dynamic constants for the airfoil. 
- T::Array{TF, 1} - A vector of floats holding the time constants for the airfoil. 
"""
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

"""
    simpleairfoil(polar)

A function that takes a simple airfoil polar to make a dynamic airfoil.
This function makes some gross assumptions.

### Inputs
- polar::Array{TF, 2} - A simple airfoil polar including alpha, lift, and drag.

### Outputs
- Airfoil

### Notes
- The function assumes that the lift curve slope is 2 pi, that the zero lift angle of attack is zero, and uses AeroDyn's default dynamic airfoil coefficients. The seperation angles of attack are the angles of attack at the minimum and maximum lift. 
"""
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


"""
    airfoil(polar; A = [0.3, 0.7], b = [0.14, 0.53], T = [1.7, 3.0])

A slightly more complex version of simpleairfoil. Takes a polar and numberically finds some characteristics. 

### Inputs
- polar::Array{TF, 2}
- A::Array{TF, 1}
- b::Array{TF, 1}
- T::Array{TF, 1}

### Outputs
- Airfoil
"""
function airfoil(polar; A = [0.3, 0.7], b = [0.14, 0.53], T = [1.7, 3.0])
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



export Functional, Iterative

abstract type DEType end #Is the model designed to be solved in one go, or iteratively, updating p ever iteration. 

struct Functional <: DEType
end

struct Iterative <: DEType
end


export Riso

abstract type DSModel end

struct Riso <: DSModel
    detype::DEType 
    n::Int #Number of airfoils simulated
    airfoils::Array{Airfoil,1}
end

function riso(airfoils; detype::DEType=Iterative()) #::Array{Airfoil,1} #TODO: I'm not sure how to type this. 
    n = length(airfoils)
    return Riso(detype, n, airfoils)
end

struct IndicialRiso <: DSModel
end

struct BeddoesLeishman <: DSModel
end

struct IndicialBeddoesLeishman <: DSModel
end