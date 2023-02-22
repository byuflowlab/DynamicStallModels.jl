using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, OpenFASTsr, LaTeXStrings

dsm = DynamicStallModels
of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)

c = 0.1

M = 0.379
a = 343.0
Vrel = M*a #60

# polar = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/learning/exploring/exampledata/NACA4412.dat", '\t'; skipstart=3) 
# af = airfoil(polar) #Todo: This constructor is broken.

du21_a17 = of.read_airfoilinput("../../data/airfoils/DU40_A17.dat") 
af = of.make_dsairfoil(du21_a17) #Todo: I think this polar might be my problem. I should try a different polar.... which means that I need to fix the constructor. :| 

# af = update_airfoil(af, A=[4.0], dcndalpha=6.320368333107256, alpha0=-0.0033903071711640564, sfun=dsm.LSP(), alphasep=af.alphasep.*5)
af = update_airfoil(af, A=[4.0], sfun=dsm.LSP(), alphasep=af.alphasep.*4.5)

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af


dsmodel = Oye(Indicial(), 1, airfoils, 1, 2)
#Note: alphasep is much higher for Faber's implemenation of the dsmodel. -> It might need more tuning... but it's something. 


tvec = range(0, 0.05, 1000) #0:0.001:0.05
Uvec = Vrel.*ones(length(tvec))

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

alphavec = alpha.(tvec)

states, loads = solve_indicial(dsmodel, [c], tvec, Uvec, alphavec)


stateplt = plot(tvec, states[:,1], leg=false, xaxis="time (s)", yaxis="f")
# display(stateplt)


cn = loads[:,1]
if dsmodel.cflag == 2
    cn_static = af.cn.(alphavec)
else
    cn_static = af.cl.(alphavec)
end

cnplt = plot( xaxis="time (s)", yaxis=L"C_n", leg=:topright)
plot!(tvec, cn, lab="DSM")
plot!(tvec, cn_static, lab="Static")
# display(cnplt) 


cyclecnplt = plot(xaxis="Angle of Attack (deg)", yaxis=L"C_n", leg=:bottomright)
plot!(alphavec.*(180/pi), cn, lab="DSM")
plot!(alphavec.*(180/pi), cn_static, lab="Static")
# vline!([af.alphasep[2]*(180/pi)], lab=L"\alpha_s")
display(cyclecnplt) 

#=
I think that I have the models matching, the difference is that I'm using Cn in this model and Cl in the other. 
=#


nothing