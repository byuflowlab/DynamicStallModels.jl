export simpleairfoil, airfoil, riso, Airfoil #TODO: I need to organize this file. 

"""
    Airfoil(polar::Array{TF, 2}, cl::Tfit, cd::Tfit, cm::Tfit, dcldalpha::TF, alpha0::TF, alphasep::Array{TF, 1}, A::Array{TF, 1}, b::Array{TF, 1}, T::Array{TF, 1})

A struct to hold all of the airfoil polars, fits, and dynamic coefficients. 

### Inputs
- polar - A matrix of floats describing the airfoil polar. This includes the angle of attack (radians), and coefficients of lift, drag, and moment. 
- cl - A fit of the coefficient of lift as a function of the angle of attack (radians). 
- cd - A fit of the coefficient of drag as a function of the angle of attack (radians). 
- cm - A fit of the coefficient of moment as a function of the angle of attack (radians). 
- dcldalpha - The lift curve slope in the linear region (1/radians). Typically near 2 pi. 
- alpha0 - The zero lift angle of attack (radians). 
- alphasep - A vector of floats holding the angles of attack at which flow fully seperates from the airfoil. In order from least to greatest. 
- A - A vector of floats holding the A dynamic constants for the airfoil. 
- b - A vector of floats holding the b dynamic constants for the airfoil. 
- T - A vector of floats holding the time constants for the airfoil. 
- S - A vector of floats holding the S constants are best fit constants for the separation point curve. 
- s - A fit of the the separation point curve. 
- xcp - The distance from the quarter chord to the center of pressure? 
"""
struct Airfoil{TF, Tfit, Fun}  
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
    sfun::Fun
    xcp::TF
end

abstract type SeparationPoint end

struct SP{fit} <:SeparationPoint #TODO: I probably don't want to call this cardoza... but whatever for now. -> Maybe just call it separation fit...
    ffit::fit
end

function SP(polar, alpha0, alphasep, dcndalpha)
    alpha = view(polar, :, 1)
    cn = view(polar, :, 2)

    fvec = reverse_separationpointcalculation(alpha, cn, dcndalpha, alpha0, alphasep)
    ffit = Akima(alpha, fvec)
    return SP(ffit)
end


struct ADSP{TF} <: SeparationPoint 
    S::Array{TF, 1}
end

struct ADFSP{fit} <: SeparationPoint
    ffit::fit
end

function ADFSP(polar, alpha0, alphasep, dcndalpha)
    alpha = view(polar, :, 1)
    cn = view(polar, :, 2)

    fvec = reverse_separationpointcalculation_ADO(alpha, cn, dcndalpha, alpha0, alphasep)
    ffit = Akima(alpha, fvec)
    return ADFSP(ffit)
end

struct BLSP{TF} <: SeparationPoint
    S::Array{TF, 1}
end

struct RSP <: SeparationPoint
end



#Note: Considering adding fit parameters.
# -> Maybe another solution would be to make an abstract type that is an Airfoil, then have structs like RisoAirfoil, and BLAirfoil. 

"""
    simpleairfoil(polar::Array{TF, 2})

A function that takes a simple airfoil polar to make a dynamic airfoil. The function assumes that the lift curve slope is 2 pi, that the zero lift angle of attack is zero, and uses AeroDyn's default dynamic airfoil coefficients. The seperation angles of attack are the angles of attack at the minimum and maximum lift.

### Inputs
- polar - A simple airfoil polar including alpha, lift, and drag.

### Outputs
- Airfoil

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

    sfun(x) = 0.0

    xcp = 0.2
    return Airfoil(polar, cl, cd, cm, dcldalpha, alpha0, alphasep, A, b, T, sfun, xcp)
end




"""
    airfoil(polar::Array{TF, 2}; A::Array{TF, 1} = [0.3, 0.7], b::Array{TF, 1} = [0.14, 0.53], T::Array{TF, 1} = [1.7, 3.0])

A slightly more complex version of simpleairfoil. Takes a polar and numerically finds some characteristics. 

### Inputs
- polar - A matrix of floats describing the airfoil polar. This includes the angle of attack (radians), and coefficients of lift, drag, and moment.
- A - A vector of floats holding the A dynamic constants for the airfoil. 
- b - A vector of floats holding the b dynamic constants for the airfoil. 
- T - A vector of floats holding the time constants for the airfoil. 
- separationpointfit - An integer telling which fit function to use. 1 -> AeroDyn Cn separation point function, 

### Outputs
- Airfoil
"""
function airfoil(polar; A = [0.3, 0.7], b = [0.14, 0.53], T = [1.7, 3.0], xcp=0.2, sfun::Union{SeparationPoint, Function}=ADFSP(1), S=zeros(4))
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

    if isa(sfun, ADFSP)
        sfun = ADFSP(polar, alpha0, alphasep, dcldalpha)
    elseif !isa(sfun, Union{Function, SeparationPoint}) #If it isn't a function or a SeparationPoint.
        @warn("The separation point function must be a function, or one of the provided options. Returning to default ADSP().")
        sfun = ADFSP(polar, alpha0, alphasep, dcldalpha)
    end

    return Airfoil(polar, cl, cd, cm, dcldalpha, alpha0, alphasep, A, b, T, sfun, xcp)
end

export update_airfoil

function update_airfoil(airfoil::Airfoil; polar=nothing, dcldalpha=nothing, alpha0=nothing, alphasep=nothing, A=nothing, b=nothing, T=nothing, sfun=nothing, xcp=nothing)

    if !(polar == nothing)
        newpolar = polar
    else
        newpolar = airfoil.polar
    end
    newcl = Akima(newpolar[:,1], newpolar[:,2])
    newcd = Akima(newpolar[:,1], newpolar[:,3])
    newcm = Akima(newpolar[:,1], newpolar[:,4])

    if !(dcldalpha==nothing)
        newslope = dcldalpha
    else
        newslope = airfoil.dcldalpha
    end

    if !(alpha0==nothing)
        newalpha0 = alpha0
    else
        newalpha0 = airfoil.alpha0
    end

    if !(alphasep==nothing)
        newalphasep = alphasep
    else
        newalphasep = airfoil.alphasep
    end

    if !(A==nothing)
        newA = A
    else
        newA = airfoil.A
    end

    if !(b==nothing)
        newb = b
    else
        newb = airfoil.b
    end

    if !(T==nothing)
        newT = T
    else
        newT = airfoil.T
    end

    if !(sfun==nothing)
        newsfun = sfun
    else
        newsfun = airfoil.sfun
    end

    if !(xcp==nothing)
        newxcp = xcp
    else
        newxcp = airfoil.xcp
    end


    return Airfoil(newpolar, newcl, newcd, newcm, newslope, newalpha0, newalphasep, newA, newb, newT, newsfun, newxcp)
end




separationpoint(airfoil::Airfoil, alpha) = separationpoint(airfoil.sfun, airfoil, alpha)

separationpoint(sfun::SP, airfoil::Airfoil, alpha) = sfun.ffit(alpha)



function reverse_separationpointcalculation(alpha, cn, dcndalpha, alpha0, alphasep)
    n = length(cn)
    fvec = Array{eltype(cn)}(undef, n)

    for i = 1:n
        sqrtarg = abs(cn[i]/(dcndalpha*(alpha[i]-alpha0)))
        f = (2*sqrt(sqrtarg) - 1)^2

        if f>1
            f=1
        elseif f<0
            f=0
        end

        if alphasep[1]<alpha[i]<alphasep[2]
            f=1
        end

        fvec[i] = f
    end

    _, alpha1idx = nearestto(alpha, alphasep[1]) # Find negative stall angle
    _, alpha2idx = nearestto(alpha, alphasep[2]) # Find positive stall angle
    
    min_positivestall, minidx_positivestall = findmin(fvec[alpha2idx:end]) # Search the positive stall region for the minimum
    min_negativestall, minidx_negativestall = findmin(fvec[1:alpha1idx]) # Search the negative stall region for the minimum

    fvec[alpha2idx+minidx_positivestall:end] .= min_positivestall # Replace all of the post positive stall region with the minimum. 
    fvec[1:minidx_negativestall] .= min_negativestall # Replace all of the post negative stall region with the minimum. 


    return fvec
end


function separationpoint(sfun::ADSP, airfoil::Airfoil, alpha)
    if isapprox(alpha, -pi, atol=1e-4)
        println("AeroDyn Original sep function called. ")
    end
    alpha0 = airfoil.alpha0
    alpha1 = airfoil.alphasep[2]
    alpha2 = airfoil.alphasep[1]
    S1, S2, S3, S4 = sfun.S

    if (alpha0 <= alpha < alpha1)
        return 1 - 0.3*exp((alpha-alpha1)/S1)
    elseif (alpha2 <= alpha < alpha0)
        return 1 - 0.3*exp((alpha2-alpha)/S3)
    elseif alpha>alpha1
        return 0.04 + 0.66*exp((alpha1-alpha)/S2)
    elseif alpha<alpha2
        return 0.04 + 0.66*exp((alpha-alpha2)/S4)
    end

    @warn("No seperation point found. Returning 1.0.")
    return 1.0
end

function reverse_separationpointcalculation_ADO(alpha, cn, dcndalpha, alpha0, alphasep)
    n = length(cn)
    fvec = Array{eltype(cn)}(undef, n)

    for i = 1:n
        sqrtarg = abs(cn[i]/(dcndalpha*(alpha[i]-alpha0)))
        f = (2*sqrt(sqrtarg) - 1)^2

        if f>1
            f=1
        elseif f<0
            f=0
        end

        if alphasep[1]<alpha[i]<alphasep[2]
            f=1
        end

        fvec[i] = f
    end
    return fvec
end

### AeroDyn separation point fit based on the static lift curve #TODO: Might need to rotate to the normal coefficient. 
separationpoint(sfun::ADFSP, airfoil::Airfoil, alpha) = sfun.ffit(alpha)

### Beddoes Leishman separation point function 
function separationpoint(sfun::BLSP, airfoil::Airfoil, alpha) 
    if isapprox(alpha, -pi, atol=1e-4)
        println("Beddoes Leishman sep function called. ")
    end

    alpha1, alpha2 = airfoil.alphasep
    S = sfun.S

    if alpha1<=alpha<=alpha2
        return 1.0-0.3*exp((alpha1-alpha)/S[1])
    elseif alpha<alpha1
        return 0.04 + 0.66*exp((alpha-alpha1)/S[2])
    else
        return 0.04 + 0.66*exp((alpha2-alpha)/S[2]) #One paper has a plus here, another has a minus. 
    end
end



# function fst_riso(alpha, liftfit, dcldalpha, alpha0)
    
#     f =  (2*sqrt(abs(liftfit(alpha)/(dcldalpha*(alpha-alpha0)))) - 1)^2
#     if f>= 1 || isnan(f)
#         return 1.0
#     elseif f<0
#         return 0.0
#     else
#         return f
#     end
    
# end

### Riso separation point function, including some corrections. 
function separationpoint(sfun::RSP, airfoil::Airfoil, alpha)
    if isapprox(alpha, -pi, atol=1e-4)
        println("Riso sep function called. ")
    end

    afm, afp = airfoil.alphasep
    clfit = airfoil.cl
    dcldalpha = airfoil.dcldalpha
    alpha0 = airfoil.alpha0
    #TODO: I'm not really sure that using the minimum of these two is really the way to avoid the problem of this blowing up to infinity. (When alpha=alpha0) (This check happens in the if statement.)

    #TODO: I'm not sure that using the absolute value here is the correct way to dodge the problem of crossing the x axis at different times.

    #Todo. I don't like using if statements for angles of attack outside of the range of attached flow. I want the function to just drive to zero. -> Todo. Plot the seperation function as a function of alpha. -> Removed the if statements and the function naturally drives to zero like the paper (not at the angles as described by Hansen.)
    
    
    ### The theory says that if the angle of attack is less (or greater than) than the seperation point aoa, then the function should return 0. -> However, using these if statements create discontinuities in the seperation point function that don't appear to be in Hansen's implementation. -> It only creates discontinuties if the incorrect values for afm and afp are used. If afm and afp are really where f(alpha)=0 then, the function won't have any discontinuities. 
    # if alpha<afm # -> TODO: If I'm getting rid of this, then I can get rid of the afm and afp arguments. 
    #     return typeof(alpha)(0)
    # elseif alpha>afp
    #     return typeof(alpha)(0)
    # end
    # error("We got here. ")
    # println(alpha)

    # @show afm, alpha, afp
    if !(afm < alpha < afp) #Check if alpha is in the bounds. 
        # println("f was set to zero. ")
        return typeof(alpha)(0)
    end
    
    cl_static = clfit(alpha)
    cl_linear = dcldalpha*(alpha-alpha0)
    f = (2*sqrt(abs(cl_static/cl_linear))-1)^2
    # println(f)  

    if f>1 #Question: What if I don't return this? I might get Inf.... or possibly NaN... but I will less likely get 1.0... which is my problem child in the seperated coefficient of lift function. -> I fixed the fully seperated coefficient of lift function... I just plugged this function inside the other and simplified. 
        return typeof(alpha)(1)
    elseif isnan(f)
        # println("f return NaN")
        return typeof(alpha)(1)
    end

    #Todo. Hansen must have some sort of switch that stops this function from reattaching when the aoa gets really high. -> like the one where you automatically set f=0 when you're outside the bounds of afm, afp
    return f
end

separationpoint(sfun::Function, airfoil::Airfoil, alpha) = sfun(alpha)



