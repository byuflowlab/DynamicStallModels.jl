#=
The SeparationPoint and Airfoil structs and convenience functions. Note that these must be defined in the same file for compilation reasons. 
=#

export simpleairfoil, airfoil, Airfoil, update_airfoil

abstract type SeparationPoint end





#AeroDyn separation point fit. 
"""
    ADSP(ffit::fit, fcfit::fit)

A struct to hold the AeroDyn separation point function. This method is a fit of the separation point function based on the separation as a function of Cl, and the chordwise forces.
"""
struct ADSP{fit} <: SeparationPoint 
    ffit::fit
    fcfit::fit
end

#Method to create an ADSP struct. 
"""
    ADSP(alpha, cn, cc, alpha0, alphasep, dcndalpha, eta)

A method to create an ADSP struct. This method is a fit of the separation point function based on the separation as a function of Cl, and the chordwise forces.
"""
function ADSP(alpha, cn, cc, alpha0, alphasep, dcndalpha, eta) #Todo: Probably needs fixing (compare to OpenFAST)

    fvec, fcvec = reverse_separationpointcalculation_ADO(alpha, cn, cc, dcndalpha, alpha0, alphasep, eta)
    ffit = Akima(alpha, fvec) 
    fcfit = Akima(alpha, fcvec)
    return ADSP(ffit, fcfit)
end


"""
    ADGSP()
AeroDyn's separation point function with Gonzalez's modifications. 
"""
struct ADGSP <: SeparationPoint #Todo: Do I need a separate struct for this?
end

#The original Beddoes-Leishman separation point function. Todo: Change the name to match Beddoes-Leishman. 
"""
    BLSP(S::Vector{TF})

A struct to hold the original AeroDyn separation point function. This function is a best fit of the separation point function, this is a direct copy of Beddoes Leishman's original implementation.

### Inputs
- S::Vector{TF} - A vector of floats holding the S constants are best fit constants for the separation point curve.
""" 
struct BLSP{TF} <: SeparationPoint 
    S::Array{TF, 1}
end



"""
    OSP()
Øye's separation point function.
"""
struct OSP <: SeparationPoint
end


"""
    RSP()
Risø's (Hansen) separation point function.
"""
struct RSP <: SeparationPoint
end




"""   
    SP(ffit::fit, fcfit::fit)
The struct to hold Cardoza2024's separation point function.

### Inputs
- ffit - A fit of the separation point function.
- fcfit - A fit of the chordwise separation point function.
"""
struct SP{fit} <:SeparationPoint 
    ffit::fit
    fcfit::fit
end

"""
    SP(alpha, cn, cc, alpha0, alphasep, dcndalpha, eta)

Cardoza2024's implementation of a separation point function. 

### Inputs
- alpha::Vector{TF} - vector of the angles of attack (show range from -180 to 180 degrees). (radians)
- cn::Vector{TF} - vector of the normal force coefficients. (Unitless)
- cc::Vector{TF} - vector of the tangential force coefficients. (Unitless)
- alpha0::TF - the zero lift angle of attack. (radians)
- alphasep::Vector{TF} - vector of the angles of attack at which flow fully separates from the airfoil. In order of negative stall, then positive stall. (radians)
- dcndalpha::TF - the normal force curve slope in the linear region. (1/radians)
- eta::TF - the separation recovery efficiency. (Unitless)
"""
function SP(alpha, cn, cc, alpha0, alphasep, dcndalpha, eta)

    fvec, fcvec = reverse_separationpointcalculation(alpha, cn, cc, dcndalpha, alpha0, alphasep, eta)
    ffit = Akima(alpha, fvec)  
    fcfit = Akima(alpha, fcvec)
    return SP(ffit, fcfit)
end






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
    model::DSModel
    polar::Array{TF, 2} #Airfoil polar - alpha (radians), cl, cd, cm
    cl::Tfit #Fit of the coefficient of lift
    cd::Tfit #Fit of the drag coefficient
    cm::Tfit #Fit of the moment coefficient
    cn::Tfit #Fit of the normal force coefficient
    cc::Tfit #Fit of the chordwise force coefficient
    dcldalpha::TF #Lift curve slope at alpha0
    dcndalpha::TF #Normal curve slope at alpha0
    alpha0::TF #The zero lift angle of attack
    alphasep #::TFV #The separation angles of attack #Todo: Figure out how to put this as a vector or an array. 
    alphacut
    cutrad::TF
    sfun::Fun #The separation point function
    c::TF #Chord length
    xcp::TF #The center of pressure
end

#################################################################
####################### Convenience Functions ###################
#################################################################

"""
    simpleairfoil(polar::Array{TF, 2})

A function that takes a simple airfoil polar to make a dynamic airfoil. The function assumes that the lift curve slope is 2 pi, that the zero lift angle of attack is zero, and uses AeroDyn's default dynamic airfoil coefficients. The seperation angles of attack are the angles of attack at the minimum and maximum lift.

### Inputs
- polar - A simple airfoil polar including alpha, lift, and drag.

### Outputs
- Airfoil

"""
function make_simpleairfoil(polar, dsmodel::DSModel, chord)
    alphavec = @view(polar[:,1])
    clvec = @view(polar[:,2])
    cdvec = @view(polar[:,3])

    cl = Akima(alphavec, clvec)
    cd = Akima(alphavec, cdvec)
    cm = Akima(alphavec, zeros(length(alphavec)))

    cnvec = @. clvec*cos(alphavec) + cdvec*sin(alphavec)
    ccvec = @. clvec*sin(alphavec) - cdvec*cos(alphavec)

    cn = Akima(alphavec, cnvec)
    cc = Akima(alphavec, ccvec)


    dcldalpha = 2*pi
    alpha0 = 0.0
    A = [0.3, 0.7, 1.0] #TODO: Maybe do tuples... Since they aren't going to change... That might make it faster. 
    b = [0.14, 0.53, 5.0]
    T = [1.7, 3.0, 0.19]

    _, minclidx = findmin(clvec) #TODO: I have the find_seperation_alpha function. 
    _, maxclidx = findmax(clvec)

    if isa(sfun, OSP) #this kind of setup is required for hermite interpolation to work properly
        alphasep = [-32*pi/180, 32*pi/180]
    else
        alphasep = [polar[minclidx, 1], polar[maxclidx,1]]
    end

    alphacut = 45*pi/180
    cutrad = 5*pi/180

    sfun(x) = 0.0

    xcp = 0.2
    return Airfoil(dsmodel, polar, cl, cd, cm, cn, cc, dcldalpha, dcldalpha, alpha0, alphasep, alphacut, cutrad, sfun, chord, xcp)
end




"""
    make_airfoil(polar::Array{TF, 2}; A::Array{TF, 1} = [0.3, 0.7], b::Array{TF, 1} = [0.14, 0.53], T::Array{TF, 1} = [1.7, 3.0])

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
function make_airfoil(polar, dsmodel::DSModel, chord; xcp=0.2, sfun::Union{SeparationPoint, Function}=ADSP(1, 1), alphacut = 45*pi/180, cutrad = 5*pi/180, eta = 1.0) 
    #TODO: Need some sort of behavior when the provided polar is too small. 

    alphavec = polar[:,1]
    clvec = polar[:,2]
    cdvec = polar[:,3]

    cl = Akima(polar[:,1], polar[:,2])
    cd = Akima(polar[:,1], polar[:,3])
    cm = Akima(polar[:,1], zeros(length(polar[:,1])))

    cnvec = @. clvec*cos(alphavec) + cdvec*sin(alphavec)
    ccvec = @. clvec*sin(alphavec) - cdvec*cos(alphavec)

    cn = Akima(alphavec, cnvec)
    cc = Akima(alphavec, ccvec) #TODO: Find dcn dalpha


    alpha0, _ = brent(cl, -0.25, 0.25)
    _, maxclidx = findmax(polar[:,2])
    _, minclidx = findmin(polar[1:maxclidx,2])

    if isa(sfun, OSP) #this kind of setup is required for hermite interpolation to work properly
        alphasep = [-32*pi/180, 32*pi/180]
    else
        alphasep = [polar[minclidx, 1], polar[maxclidx,1]]
    end

    
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

    if isa(sfun, ADSP)
        sfun = ADSP(alphavec, cnvec, ccvec, alpha0, alphasep, dcldalpha, eta)
    elseif !isa(sfun, Union{Function, SeparationPoint}) #TODO: If it isn't a function or a SeparationPoint.... Does this do anything? I enforced function and SeparationPoint in the function arguments. 
        @warn("The separation point function must be a function, or one of the provided options. Returning to default ADSP().")
        sfun = ADSP(alphavec, cnvec, ccvec, alpha0, alphasep, dcldalpha, eta)
    end

    return Airfoil(dsmodel, polar, cl, cd, cm, cn, cc, dcldalpha, dcldalpha, alpha0, alphasep, alphacut, cutrad, sfun, chord, xcp)
end



function update_airfoil(airfoil::Airfoil; dsmodel::DSModel=nothing, polar=nothing, dcldalpha=nothing, dcndalpha=nothing, alpha0=nothing, alphasep=nothing, alphacut=nothing, cutrad=nothing, sfun=nothing, chord=nothing, xcp=nothing)

    if !(dsmodel == nothing)
        newdsmodel = dsmodel
    else
        newdsmodel = airfoil.model
    end

    if !(polar == nothing)
        newpolar = polar
    else
        newpolar = airfoil.polar
    end
    newcl = Akima(newpolar[:,1], newpolar[:,2])
    newcd = Akima(newpolar[:,1], newpolar[:,3])
    newcm = Akima(newpolar[:,1], newpolar[:,4])

    cnvec = @. newpolar[:,2]*cos(newpolar[:,1]) + newpolar[:,3]*sin(newpolar[:,1])
    ccvec = @. newpolar[:,2]*sin(newpolar[:,1]) - newpolar[:,3]*cos(newpolar[:,1])

    newcn = Akima(newpolar[:,1], cnvec)
    newcc = Akima(newpolar[:,1], ccvec)

    if !(dcldalpha==nothing)
        newslope = dcldalpha
    else
        newslope = airfoil.dcldalpha
    end

    if !(dcndalpha==nothing)
        newdcndalpha = dcndalpha
    else
        newdcndalpha = airfoil.dcndalpha
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

    if !(alphacut==nothing)
        newalphacut = alphacut
    else
        newalphacut = airfoil.alphacut
    end

    if !(cutrad==nothing)
        newcutrad = cutrad
    else
        newcutrad = airfoil.cutrad
    end

    if !(sfun==nothing)
        newsfun = sfun
    else
        newsfun = airfoil.sfun
    end

    if !(chord==nothing)
        newchord = chord
    else
        newchord = airfoil.c
    end

    if !(xcp==nothing)
        newxcp = xcp
    else
        newxcp = airfoil.xcp
    end


    return Airfoil(newdsmodel, newpolar, newcl, newcd, newcm, newcn, newcc, newslope, newdcndalpha, newalpha0, newalphasep, newalphacut, newcutrad, newsfun, newchord, newxcp)
end

function Base.getproperty(obj::AbstractVector{<:Airfoil}, sym::Symbol)
    return getfield.(obj,sym)
end





###########################################################################################
##################### Separation Point functions #########################################
##########################################################################################

separationpoint(airfoil::Airfoil, alpha) = separationpoint(airfoil.sfun, airfoil, alpha)

#TODO: Can dcndalpha_circ be calculated inside the separation point function now, or does it still need to be passed in? 
separationpoint(airfoil::Airfoil, alpha, dcndalpha_circ) = separationpoint(airfoil.sfun, airfoil, alpha, dcndalpha_circ)

separationpoint(sfun::SP, airfoil::Airfoil, alpha) = sfun.ffit(alpha)

const fclimit = (1.0 + 0.2)^2

function reverse_separationpointcalculation(alpha, Cn, Cc, dcndalpha, alpha0, alphasep, eta)

    if !(length(alpha)==length(Cn)==length(Cc))
        error("reverse_separationpointcalculation(): alpha, normal coefficient, and chordwise coefficient need to all be the same length.")
    end


    n = length(Cn)
    fvec = Array{eltype(Cn)}(undef, n)
    fcvec = Array{eltype(Cn)}(undef, n)

    for i = 1:n
        sqrtarg = abs(Cn[i]/(dcndalpha*(alpha[i]-alpha0)))
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

        ### Chordwise separation point function
        bot = eta*dcndalpha*(alpha[i]-alpha0)*tan(alpha[i])
        fc = (Cc[i]/bot + 0.2)^2 #EQ 1.32b

        if fc>fclimit #OpenFAST v3.3.0 - UnsteadyAreo.f90 line 262
            fc = fclimit
        elseif isapprox(0.0, dcndalpha, atol=1e-3)
            fc = fclimit
        elseif isapprox(0.0, alpha[i], atol=1e-5)
            fc = fclimit
        elseif isapprox(alpha[i], alpha0, atol=1e-5)
            fc = fclimit
        end

        ### Force regions in unseparation region to max. 
        if alphasep[1]<alpha[i]<alphasep[2]
            fc = fclimit
        end

        fcvec[i] = fc
    end

    #Todo: I probably want to apply this same high aoa limiting on the chordwise separation point. 

    #### Limit reattachment outside of normal angles of attack. 
    _, alpha1idx = nearestto(alpha, alphasep[1]) # Find negative stall angle
    _, alpha2idx = nearestto(alpha, alphasep[2]) # Find positive stall angle
    
    min_positivestall, minidx_positivestall = findmin(fvec[alpha2idx:end]) # Search the positive stall region for the minimum
    min_negativestall, minidx_negativestall = findmin(fvec[1:alpha1idx]) # Search the negative stall region for the minimum

    fvec[alpha2idx+minidx_positivestall:end] .= min_positivestall # Replace all of the post positive stall region with the minimum. 
    fvec[1:minidx_negativestall] .= min_negativestall # Replace all of the post negative stall region with the minimum. 


    ##### Same as above but for chordwise separation point function. 
    min_positivestall, minidx_positivestall = findmin(fcvec[alpha2idx:end]) # Search the positive stall region for the minimum
    min_negativestall, minidx_negativestall = findmin(fcvec[1:alpha1idx]) # Search the negative stall region for the minimum

    fcvec[alpha2idx+minidx_positivestall:end] .= min_positivestall # Replace all of the post positive stall region with the minimum. 
    fcvec[1:minidx_negativestall] .= min_negativestall # Replace all of the post negative stall region with the minimum. 


    


    return fvec, fcvec
end


function separationpoint(sfun::BLSP, airfoil::Airfoil, alpha)
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

# ### Beddoes Leishman separation point function 
# function separationpoint(sfun::BLSP, airfoil::Airfoil, alpha) 
#     if isapprox(alpha, -pi, atol=1e-4)
#         println("Beddoes Leishman sep function called. ")
#     end

#     alpha1, alpha2 = airfoil.alphasep
#     S = sfun.S

#     if alpha1<=alpha<=alpha2
#         return 1.0-0.3*exp((alpha1-alpha)/S[1])
#     elseif alpha<alpha1
#         return 0.04 + 0.66*exp((alpha-alpha1)/S[2])
#     else
#         return 0.04 + 0.66*exp((alpha2-alpha)/S[2]) #One paper has a plus here, another has a minus. 
#     end
# end

function reverse_separationpointcalculation_ADO(alpha, Cn, Cc, dcndalpha, alpha0, alphasep, eta)

    n = length(Cn)
    fvec = Array{eltype(Cn)}(undef, n)
    fcvec = Array{eltype(Cn)}(undef, n)

    for i = 1:n
        sqrtarg = abs(Cn[i]/(dcndalpha*(alpha[i]-alpha0)))
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

        ### Chordwise separation point function
        # bot = eta*dcndalpha*(alpha[i]-alpha0)*tan(alpha[i])
        bot = eta*dcndalpha*(alpha[i]-alpha0)*alpha[i]
        fc = (Cc[i]/bot + 0.2)^2 #EQ 1.32b

        if fc>fclimit #OpenFAST v3.3.0 - UnsteadyAreo.f90 line 262
            fc = fclimit
        end

        fcvec[i] = fc
    end

    return fvec, fcvec
end

### AeroDyn separation point fit based on the static lift curve #TODO: Might need to rotate to the normal coefficient. 
separationpoint(sfun::ADSP, airfoil::Airfoil, alpha) = sfun.ffit(alpha) #Todo: These other separation point functions need the chordwise separation point function. 

function separationpoint(sfun::ADGSP, airfoil::Airfoil, alpha)

    Cn = airfoil.cl(alpha)*cos(alpha) + (airfoil.cd(alpha) - airfoil.cd(airfoil.alpha0))*sin(alpha)

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



    afm, _ = airfoil.alphasep
    clfit = airfoil.cl
    dcldalpha = airfoil.dcldalpha
    alpha0 = airfoil.alpha0

    f(x) = clfit(x) - dcldalpha*(x - alpha0)/4

    afp , _ = brent(f , 0.17 , 0.8726)


    #println(afp)
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
        return typeof(alpha)(1.0)
    elseif isnan(f)
        # println("f return NaN")
        return typeof(alpha)(1.0)
    end
    


    #Todo. Hansen must have some sort of switch that stops this function from reattaching when the aoa gets really high. -> like the one where you automatically set f=0 when you're outside the bounds of afm, afp
    return f
end

separationpoint(sfun::Function, airfoil::Airfoil, alpha) = sfun(alpha)

#=
Larsen's separation point function from his 2007 paper. 
=#
function separationpoint(sfun::OSP, airfoil::Airfoil, alpha)
    #println("using the OSP separation point function. Currently at line 639 in airfoils.jl ")
    if alpha>airfoil.alphasep[2] #? right after stall is fully separated? not partially?
        return 0.0
    elseif alpha<airfoil.alpha0
        
        return 0.0
    else
        alpha0 = airfoil.alpha0
        cn = airfoil.cl(alpha) #Todo: I need to make this be able to switch between cl and cn. 

        cn_inv = airfoil.dcldalpha*(alpha-alpha0)
        cn_fs = cl_fullysep_faber(airfoil, alpha)

        fst = (cn - cn_fs)/(cn_inv - cn_fs)

        #define tolerances
        tol1 = .35
        tol2 = .35

        if (cn - cn_fs) < tol1 && (cn_inv - cn_fs) < tol2
            #println("fst: ", fst, ", top: (cn-cn_fs): ", (cn-cn_fs), ", bottom: (cn_inv-cn_fs): ", (cn_inv-cn_fs), ", alpha ", alpha *180/pi)
            return 1.0
        end

        if fst>1
            return 1.0  
        elseif fst<0
            return #0.0
        else
            return fst #! is something up with these values? running a for loop with alphavec on oyecomparer.jl doesn't go above .02ish
        end
    end
end

#=
Gonzalez modification of the separation point function, as found in OpenFAST v3.3.0
=#
function separationpoint(sfun::ADGSP, airfoil::Airfoil, alpha, dcndalpha_circ)
    delalpha = alpha-airfoil.alpha0

    # Cd0 = airfoil.cd(airfoil.alpha0)
    Cd0 = airfoil.model.Cd0
    # @show alpha
    
    Cn = airfoil.cl(alpha)*cos(alpha) + (airfoil.cd(alpha) - Cd0)*sin(alpha)

    if isapprox(dcndalpha_circ, 0.0)
        tr = 0
    elseif alpha == airfoil.alpha0
        tr = 0
    elseif Cn == 0
        tr = 0
    else
        tr = Cn/(dcndalpha_circ*delalpha)
        if tr<0
            tr=0
        end
    end

    fst = ((3*sqrt(tr) - 1)/2)^2

    if fst>1
        fst = 1
    end

    return fst
end






chordwiseseparationpoint(airfoil::Airfoil, alpha) = chordwiseseparationpoint(airfoil.sfun, airfoil, alpha)

chordwiseseparationpoint(airfoil::Airfoil, alpha, dcndalpha_circ) = chordwiseseparationpoint(airfoil.sfun, airfoil, alpha, dcndalpha_circ)

chordwiseseparationpoint(sfun::SP, airfoil::Airfoil, alpha) = sfun.fcfit(alpha)

chordwiseseparationpoint(sfun::ADSP, airfoil::Airfoil, alpha) = sfun.fcfit(alpha)

function chordwiseseparationpoint(sfun::BLSP, airfoil::Airfoil, alpha)
    return separationpoint(airfoil, alpha)
end

function chordwiseseparationpoint(sfun::OSP, airfoil::Airfoil, alpha)
    return separationpoint(airfoil, alpha)
end

function chordwiseseparationpoint(sfun::RSP, airfoil::Airfoil, alpha)
    return separationpoint(airfoil, alpha)
end

#=
Gonzalez modification of the separation point function, as found in OpenFAST v3.3.0
=#
function chordwiseseparationpoint(sfun::ADGSP, airfoil::Airfoil, alpha, dcndalpha_circ)
    # @show alpha, dcndalpha_circ
    Cd0 = airfoil.model.Cd0

    Cc = airfoil.cl(alpha)*sin(alpha) - (airfoil.cd(alpha) - Cd0)*cos(alpha)

    # println(airfoil.cl(alpha), ",  ", airfoil.cd(alpha), ", ", Cd0, ", ", alpha)

    D = airfoil.model.eta*dcndalpha_circ*(alpha-airfoil.alpha0)*alpha

    fc = (Cc/D + 0.2)^2 #Todo. This is always over 1.44

    # delalpha = airfoil.alpha0-alpha #Note: I don't have the add_sub_2pi() function here (to bring the difference of the angles within a difference of pi), but... It doesn't look like that's the current problem. It looks like it is something else. 
    # @show delalpha

    # @show fc
    # @show Cc, D #Todo. These values are off, not by tons, but off. -> Cd0 was off. 

    return min(fc, fclimit)
end