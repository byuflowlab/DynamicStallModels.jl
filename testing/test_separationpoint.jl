using DynamicStallModels, OpenFASTsr, FLOWMath, Plots, Plots.PlotMeasures, DelimitedFiles, LaTeXStrings

#=
Test the separation point functions. (Starting with the AeroDyn original separation point fit (Equation 1.32 inversed.))

Adam present 9/27/22

Chordwise separation point functions added. - Adam 11/28/22
=#

of = OpenFASTsr
dsm = DynamicStallModels

path = dirname(@__FILE__)
cd(path)




### Create airfoil
c = 1.0

du21_a17 = of.read_airfoilinput("../data/DU21_A17.dat")


af = of.make_dsairfoil(du21_a17)


# af = update_airfoil(af; dcldalpha=af.dcldalpha*1.2) 



alpha = (-180:1:180).*(pi/180)
clvec = af.cl(alpha)
cdvec = af.cd(alpha)
cnvec = af.cn(alpha)
ccvec = af.cc(alpha)
aoavec = alpha.*(180/pi)

linearlift(airfoil, aoa) = airfoil.dcldalpha*(aoa-airfoil.alpha0)
lcl = linearlift.(Ref(af), alpha)

farclplt = plot(aoavec, af.cl(alpha), xaxis="Angle of Attack (deg)", yaxis="Coefficient of Lift", lab="Static Data", ylims=(-5,5))
plot!(aoavec, lcl, lab="Linear data")

closeclplt = plot(aoavec, af.cl(alpha), xaxis="Angle of Attack (deg)", yaxis="Coefficient of Lift", lab="Static Data", xlims=(-15,15), ylims=(-5, 5), legend=:topleft)
plot!(aoavec, lcl, lab="Linear data")
vline!([af.alphasep].*(180/pi), lab="Alpha sep")

clplt = plot(farclplt, closeclplt, layout=(2,1))
display(clplt)





#TODO: Change the name to just SP, so it doesn't have my name. 
#Todo: I wonder why there is a difference betewen the ADO fit and the my fit. 
sp = dsm.SP(alpha, cnvec, ccvec, af.alpha0, af.alphasep, af.dcldalpha, af.eta)
f_present = dsm.separationpoint.(Ref(sp), Ref(af), alpha)
fc_present = dsm.chordwiseseparationpoint.(Ref(sp), Ref(af), alpha)

# S = [du21_a17.s1, du21_a17.s2, du21_a17.s3, du21_a17.s4]
# ado = dsm.BLSP(S) 
# f_ADO = dsm.separationpoint.(Ref(ado), Ref(af), alpha)  


# ffitvec_ADO = dsm.separationpoint.(Ref(af.sfun), Ref(af), alpha)
# fcfitvec_ADO = dsm.chordwiseseparationpoint.(Ref(af.sfun), Ref(af), alpha)

S = [du21_a17.s1, du21_a17.s2, du21_a17.s3, du21_a17.s4]
BLsep = dsm.BLSP(S)#S coefficients from NACA 0012
f_BL = dsm.separationpoint.(Ref(BLsep), Ref(af), alpha) 

risosep = dsm.RSP()
f_riso = dsm.separationpoint.(Ref(risosep), Ref(af), alpha) 


alphasep_riso = sort([dsm.find_seperation_alpha(af.cl, af.dcldalpha, af.alpha0)...]) #Todo: I need to consider having this as the default for the riso. 
af_riso = update_airfoil(af; alphasep=alphasep_riso)
f_riso_risosep = dsm.separationpoint.(Ref(risosep), Ref(af_riso), alpha) 

larsensep = dsm.LSP()
f_larsen = dsm.separationpoint.(Ref(larsensep), Ref(af), alpha)


f1plt = plot( xaxis="Angle of Attack (deg)", yaxis="Separation Point f", leg=:topleft, xlims=(-180,10))
# plot!(aoavec, ffitvec_ADO, lab="fit ADO")
plot!(aoavec, f_present, lab="fit present")
# plot!(aoavec, f_ADO, lab="f ADO")
plot!(aoavec, f_BL, lab="BL")
plot!(aoavec, f_riso, lab="Risø, ADO alphasep", markershape=:cross)
plot!(aoavec, f_riso_risosep, lab="Risø, Riso alphasep", markershape=:x)
plot!(aoavec, f_larsen, lab="Larsen")
vline!([af.alphasep].*(180/pi), lab="Alpha sep")

f2plt = plot( xaxis="Angle of Attack (deg)", yaxis="Separation Point f", leg=false, xlims=(-20,180))
# plot!(aoavec, ffitvec_ADO, lab="fit ADO")
# plot!(aoavec, f_present, lab="fit present")
# plot!(aoavec, f_ADO, lab="f ADO")
# plot!(aoavec, f_BL, lab="f BL")
# plot!(aoavec, f_riso, lab="f_riso, ADO alphasep",markershape=:cross)
# plot!(aoavec, f_riso_risosep, lab="f_riso, Riso alphasep", markershape=:x)
plot!(aoavec, f_larsen, lab="Larsen")
vline!([af.alphasep].*(180/pi), lab="Alpha sep")
vline!([af.alpha0*180/pi], lab="alpha0")

fplt = plot(f1plt, f2plt, layout=(2,1))
display(fplt)
# savefig("fsepplt.png")

#=
f_ADO is so jacked up because the S constants are all zero. 

Well f_BL is jacked up as well. It's a little more jacked than I expected. Actually, it's not too terrible, but it doesn't ever reattach. And it needs to be augmented to drop off on the other side. 

Something must be broken with sparationpoint_riso. How odd... even using the alphasep points it's still bad... Honestly, I'm a little confused why though... because fst is a simpler version of the same equation... So why is it garbage? 


#Note: I think how I want to implement this is precalculate all of the f values for all of the possible alphas (create a fit), then just use the fit at calc time. So when I create the airfoil, I decide what kind of separation point fit function I want to use. 

I think that I realized why fit ADO and fit present are different. The difference is the input dcldalpha. I update the dcldalpha, but the ADO function gets initialized with the old dcldalpha. 

=#

f3plt = plot( xaxis="Angle of Attack (deg)", yaxis=L"f'_c", leg=:topleft, xlims=(-180,10))
# plot!(aoavec, fcfitvec_ADO, lab="fit ADO")
plot!(aoavec, fc_present, lab="fit present")
# plot!(aoavec, f_ADO, lab="f ADO")
# plot!(aoavec, f_BL, lab="f BL")
# plot!(aoavec, f_riso, lab="f_riso, ADO alphasep", markershape=:cross)
# plot!(aoavec, f_riso_risosep, lab="f_riso, Riso alphasep", markershape=:x)
vline!([af.alphasep].*(180/pi), lab="Alpha sep")

f4plt = plot( xaxis="Angle of Attack (deg)", yaxis=L"f'_c", leg=false, xlims=(-20,180))
# plot!(aoavec, fcfitvec_ADO, lab="fit ADO")
plot!(aoavec, fc_present, lab="fit present")
# plot!(aoavec, f_ADO, lab="f ADO")
# plot!(aoavec, f_BL, lab="f BL")
# plot!(aoavec, f_riso, lab="f_riso, ADO alphasep",markershape=:cross)
# plot!(aoavec, f_riso_risosep, lab="f_riso, Riso alphasep", markershape=:x)
vline!([af.alphasep].*(180/pi), lab="Alpha sep")

fcplt = plot(f3plt, f4plt, layout=(2,1))
# display(fcplt)

#=
 I was redefining dcndalpha before after I initiated the fit function for the ADO, so my fit was different. If I don't do that, then the chordwise separation point functions match (for now). 

I don't expect the little dip in the center. I wonder if AeroDyn does that as well. I'll just have to compare. 
=#







nothing