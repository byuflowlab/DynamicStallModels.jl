using DelimitedFiles, Plots, FLOWMath, Polynomials, Roots, Dierckx
pnoms = Polynomials

function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

function cap!(xvec, capval)
    for i=1:length(xvec)
        if xvec[i]>capval
            xvec[i] = capval
        end
    end
    return xvec
end

function button!(xvec, buttval)
    for i=1:length(xvec)
        if xvec[i]<buttval
            xvec[i] = buttval
        end
    end
    return xvec
end

polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989staticdata_fig1.csv", ',')

polar[:,1] = polar[:,1].*(pi/180)

Clfit = Akima(polar[:,1], polar[:,2])

extpolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/extendedNACA0012_fromexperimental.dat", ',')
extClfit = Akima(extpolar[:,1], extpolar[:,2])

Clpolar = scatter(polar[:,1].*(180/pi), polar[:,2], label="Original")
plot!(extpolar[:,1].*(180/pi), extpolar[:,2], label="Extended")

zeroalf, zeroidx = nearestto(polar[:,2], 0.0) #Nearest to zero lift
tenalf, tenidx = nearestto(polar[:,1], 10.0*pi/180) #Ten is typically   before stall
linrng = zeroidx:tenidx
mxb = pnoms.fit(polar[linrng,1], polar[linrng,2], 1)
dcldalpha = mxb.coeffs[2]#least squares fit
alpha0 = pnoms.roots(mxb)[1] #least squares fit

### using the extended data has shown how what I was doing before wouldn't always get what I wanted.
## The zero lift should only occur once in the range (-50, 50)
alfn50, alfn50idx = nearestto(extpolar[:,1], -50*pi/180) #Relative to the entire vector
alf50, alf50idx = nearestto(extpolar[:,1], 50*pi/180) #Relative to the entire vector
extzeroCl, extzeroidx = nearestto(extpolar[alfn50idx:alf50idx,2], 0.0) #Nearest to zero lift #Relative to a subset
extzeroidx = extzeroidx + alfn50idx - 1  #Should be relative to the whole vector
extzeroalf = extpolar[extzeroidx, 1] #Todo: I'm not sure I need this line. 

extmaxcl, extmaxclidx = findmax(extpolar[extzeroidx:alf50idx,2]) #Relative to a subset
extmaxclidx = extmaxclidx + extzeroidx - 1 #Should be relative to the whole vector
extmaxcl5, extmaxclidx5 = nearestto(extpolar[extzeroidx:extmaxclidx,2], 0.5*extmaxcl) #Half of the Cl of max Cl
extmaxclidx5 = extmaxclidx5 + extzeroidx - 1 #Should be relative to the whole vector
extlinrng = extzeroidx:extmaxclidx5
extmxb = pnoms.fit(extpolar[extlinrng,1], extpolar[extlinrng,2], 1)
extdcldalpha = extmxb.coeffs[2]#least squares fit
extalpha0 = pnoms.roots(extmxb)[1] #least squares fit
extalpha_maxcl = extpolar[extmaxclidx, 1]

# println("-50 deg: ", alfn50*180/pi, "   idx: ", alfn50idx)
# println(" 50 deg: ", alf50*180/pi, "   idx: ", alf50idx)
# println(" Cl0: ", extzeroCl, "   idx: ", extzeroidx, "   val: ", extpolar[extzeroidx, 2])
# println(" max Cl: ", extmaxcl, "   idx: ", extmaxclidx, "   val: ", extpolar[extmaxclidx, 2])
# println(" 1/2 Cl: ", extmaxcl5, "   idx: ", extmaxclidx5, "   val: ", extpolar[extmaxclidx5, 2], "  anlytical: ", extmaxcl/2)
# println("  range: ", extlinrng)

println("")
println("dcldalpha: ", dcldalpha)
println("alpha0: ", alpha0)

println("")
println("Extended data:")
println("  max Cl: ", extmaxcl)
println("  dcldcalpha: ", extdcldalpha)
println("  alpha0: ", extalpha0)

tempalfvec = collect(-6:0.05:14).*(pi/180)
Clf = extmxb.(tempalfvec)

Clpolar = scatter(polar[:,1].*(180/pi), polar[:,2], label="Original", legend=:topleft)
plot!(extpolar[:,1].*(180/pi), extpolar[:,2], label="Extended")
plot!(extpolar[extlinrng,1].*(180/pi), extpolar[extlinrng,2], label="rng")
plot!(tempalfvec.*(180/pi), Clf, label="fit")
xlims!((-6, 14))
display(Clpolar)

function residual(alpha)
    Cl = Clfit(alpha)
    return (2*sqrt(Cl/(dcldalpha*(alpha-alpha0)))-1)^2 - 0.7
end

alpha7 = find_zero(residual, 10*pi/180) 

function extresidual(alpha)
    Cl = extClfit(alpha)
    return (2*sqrt(Cl/(extdcldalpha*(alpha-extalpha0)))-1)^2 - 0.7
end
extalpha7 = find_zero(extresidual, (9*pi/180, 20*pi/180), Bisection()) 

nothing

function residual2(alfa)
    Cl = Clfit(alfa)
    return (dcldalpha*(alfa-alpha0)/4)-Cl
end

alfvec = collect(alpha7:0.0001:30*pi/180)
r2vec = residual2.(alfvec)

function extresidual2(alfa)
    Cl = extClfit(alfa)
    return (extdcldalpha*(alfa-extalpha0)/4)-Cl
end

extalphasep = find_zero(extresidual2, (extalpha7, 5*extalpha7), Bisection())

extalfvec = collect(extalpha7:0.0001:30*pi/180)
extr2vec = extresidual2.(extalfvec)

plt = plot(alfvec, r2vec, label="original")
plot!(extalfvec, extr2vec, label="extended")
xlabel!("Angle of attack (radians)")
ylabel!("Residual function value")

display(plt)

extinviscidlift(alf) = extdcldalpha*(alf-extalpha0)

delCl0 = abs(extinviscidlift(extalpha0)-extClfit(extalpha0))
delClstall = abs(extinviscidlift(extalpha7)-extClfit(extalpha7)) #extalpha_maxcl
m = -0.3/(delClstall-delCl0)
fClfit_linear(delCl) = m*delCl + 1.0



delClsep = abs(extinviscidlift(extalphasep)-extClfit(extalphasep))
fClfit_Akima = Akima([delCl0, delClstall, delClsep], [1.0, 0.7, 0.0])
fClfit_poly2 = Spline1D([delCl0, delClstall, delClsep], [1.0, 0.7, 0.0], k=2, s=0.0)
m2 = -1.0/(delClsep-delClstall)
fClfit_linear2(delCl) = m2*delCl + 1.0

function delCl(alpha)
    Cl = extClfit(alpha)
    Cli = extinviscidlift(alpha)
    return abs(Cli-Cl)
end

function ffunk_linear(alpha)
    Cl = extClfit(alpha)
    Cli = extinviscidlift(alpha)
    return fClfit_linear(abs(Cli-Cl))
end

function ffunk_linear2(alpha)
    Cl = extClfit(alpha)
    Cli = extinviscidlift(alpha)
    f = fClfit_linear2(abs(Cli-Cl))
    if f>1.0
        f = 1.0
    elseif f<0.0
        f = 0.0
    end
    return f
end

function ffunk_linear3(alpha; sig=false)
    #This causes a large discontinuity. It needs a smoothing function between the two functions. Which means the whole value needs to be calculated... then splined. :(
    Cl = extClfit(alpha)
    Cli = extinviscidlift(alpha)
    
    if alpha<extalpha7
        f = fClfit_linear(abs(Cli-Cl))
    else
        f = fClfit_linear2(abs(Cli-Cl))
    end  
    
    if f>1.0
        f = 1.0
    elseif f<0.0
        f = 0.0
    end
    return f
end

function ffunk_Akima(alpha)
    Cl = extClfit(alpha)
    Cli = extinviscidlift(alpha)
    return fClfit_Akima(abs(Cli-Cl))
end

function ffunk_poly2(alpha)
    Cl = extClfit(alpha)
    Cli = extinviscidlift(alpha)
    return fClfit_poly2(abs(Cli-Cl))
end

function hansen(alpha)
    Cl = Clfit(alpha)
    return (2*sqrt(Cl/(dcldalpha*(alpha-alpha0)))-1)^2 
end

function hansen_ext(alpha)
    Cl = extClfit(alpha)
    f = (2*sqrt(Cl/(extdcldalpha*(alpha-extalpha0)))-1)^2
    if f>1.0
        f = 1.0
    end
    return f
end

function leishman(alphan, alpha1, alpha0) #We could try some other f functions.
    S = [0.025, 0.025] #I got this from trying to match the static data 
    if alphan<=alpha1
        return 1.0 - 0.3*exp((alphan-alpha1)/S[1])
    else
        return 0.04 + 0.66*exp((alpha1-alphan)/S[2])
    end
end





alfvec = collect(0.0:0.0001:40*pi/180)
aoavec = alfvec.*(180/pi)
hvec = hansen.(alfvec)
hextvec = hansen_ext.(alfvec)
lvec = leishman.(alfvec, Ref(alpha7), Ref(alpha0))
flinearvec = ffunk_linear.(alfvec)
flinear2vec = ffunk_linear2.(alfvec)
flinear3vec = ffunk_linear3.(alfvec)
# fsig = ffunk_linear3.(alfvec;sig=true)
fakimavec = ffunk_Akima.(alfvec)
fpoly2vec = ffunk_poly2.(alfvec)

delClvec = delCl.(alfvec)
fCllin = fClfit_linear.(delClvec)
fCllin2 = fClfit_linear2.(delClvec)
delClt = delCl(extalpha7)
hardness = 4
fsig = sigmoid_blend.(fCllin, fCllin2, delClvec, delClt, hardness)
fsig = button!(fsig, 0.0)



plt2 = plot(aoavec, hvec, legend=:outerright, label="Hansen", size=(1250,600))
plot!(aoavec, hextvec, label="Hansen - Extended Data")
plot!(aoavec, lvec, label="Leishman")
plot!(aoavec, flinearvec, label="f fit - Linear")
plot!(aoavec, flinear2vec, label="f fit - Linear 2", linestyle=:dash)
plot!(aoavec, flinear3vec, label="f fit - Linear 3", linestyle=:dot)
plot!(aoavec, fsig; label="f fit - Linear 3 Sigmoid", linestyle=:dashdotdot)
plot!(aoavec, fakimavec, label="f fit - Akima")
plot!(aoavec, fpoly2vec, label="f fit - poly 2")
vline!([alpha7*180/pi], label="α7", linealpha=0.3)
vline!([alpha0*180/pi], label="α0", linealpha=0.3)
vline!([extalphasep*180/pi], label="αsep", linealpha=0.3)
ylabel!("seperation point")
xlabel!("angle of attack (degrees)")
ylims!((-0.5, 1.25))
display(plt2)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/beddoesleishmanlarsen/separation/separationalgs.png")

# alpha7 = find_zero(residual3, 10*pi/180) 

inviscidlift(alf) = dcldalpha*(alf-alpha0)
extinviscidlift(alf) = extdcldalpha*(alf-extalpha0)

alfvec2 = collect(0.0:0.0001:35*pi/180)
aoavec2 = alfvec2.*(180/pi)

cli = inviscidlift.(alfvec2)
clsf = Clfit.(alfvec2)

extcli = extinviscidlift.(alfvec2)
extclsf = extClfit.(alfvec2) 

staticpolar = plot(polar[:,1].*(180/pi), polar[:,2], legend=:topleft, label="Static")
plot!(aoavec2, clsf, label="Static fit")
plot!(aoavec2, extclsf, label="Static fit - Extended ")
plot!(aoavec2, cli, label="Inviscid")
plot!(aoavec2, extcli, label="Invscid - Extended")
plot!(aoavec2, cli./4, label="Inviscid/4")
plot!(aoavec2, extcli./4, label="Inviscid/4 Extended")
vline!([extalphasep*180/pi], label="αsep")
ylims!((-1.0, 7.5))
xlabel!("Angle of attack (degrees)")
ylabel!("Coefficient of Lift")
display(staticpolar)


