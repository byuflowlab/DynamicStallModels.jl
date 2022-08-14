using Plots
using DelimitedFiles
using FLOWMath
using Polynomials

include("Oye.jl")

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

polar = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989staticdata_fig1.csv", ',')
polar[:,1] = polar[:,1].*(pi/180)
cl = Akima(polar[:,1], polar[:,2])

function alpha(t)
    c = 0.1
    M = 0.379
    a = 343.0
    shift = 10.3
    amp = 8.1
    k = 0.075

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

function U(t)
    M = 0.379
    a = 343.0
    return M*a
end

zeroalf, zeroidx = nearestto(polar[:,2], 0.0) #Nearest to zero lift
tenalf, tenidx = nearestto(polar[:,1], 10.0*pi/180) #Ten is typically before stall #*pi/180
linrng = zeroidx:tenidx

mxb = fit(polar[linrng,1], polar[linrng,2], 1)

dcldalpha = mxb.coeffs[2] #least squares fit 

alpha0 = roots(mxb)[1] #least squares fit   # NACA 0012 is a symmetric airfoil... so shouldn't the zero lift be at zero? 

# polarplt = plot(polar[:,1], polar[:,2], leg=false)
# plot!(mxb)
# display(polarplt)

maxcl, maxclidx = findmax(polar[:,2]) #critical lift coefficient (given by Larsen 2007), the beginning of LE separation. 
# println("Static Max Cl: ", maxcl)
alphav = polar[maxclidx, 1] 

maxcl = dcldalpha*(alphav-alpha0) #Projecting the inviscid lift to the static stall point
maxcl = 1.2

function clinv(alpha) #NOTE: I'm fudging the slope at zero lift value to match the experimental rather than calculating it
    return dcldalpha*(alpha-alpha0)
end

fs = 1.0 #fst #Assuming that we are starting at the static separation f, but we'll play
c = 0.1
M = 0.379
a = 343.0
Vrel = M*a

tspan = (0.0, 0.05)
c = 0.1

Cld, tvec, fvec = Oye(tspan, alpha, U, c; nt=100)

alphavec = alpha.(tvec)
alfavec = alphavec.*(180/pi)
Clst = cl.(alphavec)

mat = hcat(alfavec, Cld)
# writedlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/outputs/leishman1989/naca0012_8.1_10.3s_Cdn_Oye.csv", mat, ',')

exppolar = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989fig8Cn.csv", ',')


cycleplt = plot(leg=:topleft)
plot!(alfavec, Cld, lab="Oye")
plot!(alfavec, Clst, lab="Static", markershape=:cross)
scatter!(exppolar[:,1], exppolar[:,2], lab="Experimental")
display(cycleplt)

timeplt = plot(leg=:bottomleft)
plot!(tvec, Cld, lab="Oye")
plot!(tvec, Clst, lab="Static")
# display(timeplt)

fplt = plot(tvec, fvec, leg=false, xaxis="time (seconds)", yaxis="Separation - f")
# display(fplt)

nothing