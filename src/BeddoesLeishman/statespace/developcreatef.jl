using DelimitedFiles, Plots, FLOWMath, Polynomials, Roots
pnoms = Polynomials

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989staticdata_fig1.csv", ',')

polar[:,1] = polar[:,1].*(pi/180)

Clfit = Akima(polar[:,1], polar[:,2])

zeroalf, zeroidx = nearestto(polar[:,2], 0.0) #Nearest to zero lift
tenalf, tenidx = nearestto(polar[:,1], 10.0*pi/180) #Ten is typically   before stall
linrng = zeroidx:tenidx
mxb = pnoms.fit(polar[linrng,1], polar[linrng,2], 1)
dcldalpha = mxb.coeffs[2]#least squares fit
alpha0 = pnoms.roots(mxb)[1] #least squares fit

function residual(alpha)
    Cl = Clfit(alpha)
    return (2*sqrt(Cl/(dcldalpha*(alpha-alpha0)))-1)^2 - 0.7
end

alpha7 = find_zero(residual, 10*pi/180) #Note: This can be replaced by just finding the static stall angle. I don't know which one is better, but I'm pretty sure that finding the static stall angle is more reliable. 

#Assuming that f=1.0 at alpha0, f=0.7 @ alpha7. I need a third point to pin to. Hansen 2004 says that's when C_L_static = dCldalpha*(alpha-alpha0)

# I just realized, we'll also want to do something for the negative side. 

function residual2(alfa)
    Cl = Clfit(alfa)
    return (dcldalpha*(alfa-alpha0)/4)-Cl
end

alfvec = collect(alpha7:0.0001:30*pi/180)
r2vec = residual2.(alfvec)

plt = plot(alfvec, r2vec)
# display(plt)

#Well


# alphasep = find_zero(residual2, alpha7)

# function ffunc(alfa, alfa1, alfa0)
#     # println(alfa)
#     Cl = Clfit(alfa)
#     # println(Cl)
#     f = (2*sqrt(abs(Cl/(dcldalpha*(alfa-alfa0))))-1)^2
#     println(f)
#     return f
# end #Hansen 2004 pg 13 EQ 15

nothing

# function createf(polar)

#     #Todo: This function only determines the static stall seperation point... but not the fit. 
#     # I have a couple options of what I can do:
#     #1- I can assume that alpha0 correlates with f=1 and alpha1 correlates with f=0.7, then let that trend continue.... or something like that
#     #2- Another option is to use Hansen's equation to solve for f at every given alpha. Problem is that f would shoot off infinity near alpha0. Note that Hansen uses the maximum dcldalpha. 
#     # I could also just solve for alpha when f=0.7.... I'm not sure. 
#     ### Todo: I need to think about this. And read. And write it out on the whiteboard. 
#     if maximum(polar[:,1])>2*pi
#         warning("BLL: createf function requires the polar in radians.")
#     end

#     Clfit = Akima(polar[:,1], polar[:,2])

#     zeroalf, zeroidx = nearestto(polar[:,2], 0.0) #Nearest to zero lift
#     tenalf, tenidx = nearestto(polar[:,1], 10.0*pi/180) #Ten is typically   before stall
#     linrng = zeroidx:tenidx
#     mxb = fit(polar[linrng,1], polar[linrng,2], 1)
#     dcldalpha = mxb.coeffs[2]#least squares fit
#     alpha0 = roots(mxb)[1] #least squares fit

#     function ffunc(alfa, alfa1, alfa0)
#         # println(alfa)
#         Cl = Clfit(alfa)
#         # println(Cl)
#         f = (2*sqrt(abs(Cl/(dcldalpha*(alfa-alfa0))))-1)^2
#         println(f)
#         return f
#     end #Hansen 2004 pg 13 EQ 15

#     return ffunc
# end