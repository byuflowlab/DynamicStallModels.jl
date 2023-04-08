using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, OpenFASTsr, LaTeXStrings

dsm = DynamicStallModels
of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)

#original c = 0.55, M = 0.07

c = 0.55

M = 0.11
a = 343.0
Vrel = M*a #60

# polar = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/learning/exploring/exampledata/NACA4412.dat", '\t'; skipstart=3) 
# af = airfoil(polar) #Todo: This constructor is broken.

#du21_a17 = of.read_airfoilinput("../../data/airfoils/DU40_A17.dat") 
#af = of.make_dsairfoil(du21_a17) #Todo: I think this polar might be my problem. I should try a different polar.... which means that I need to fix the constructor. :| 

# af = update_airfoil(af, A=[4.0], dcndalpha=6.320368333107256, alpha0=-0.0033903071711640564)

af = airfoil(polar_0015; A=[8], sfun=dsm.LSP()) 



airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af


dsmodel = Oye(Indicial(), 1, airfoils, 1, 2)


tvec = range(0, 1.0, 1000) #0:0.001:0.05
Uvec = Vrel.*ones(length(tvec))

function alpha(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.0
    k = 0.051

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
cn_static = af.cn.(alphavec)

cnplt = plot( xaxis="time (s)", yaxis=L"C_l", leg=:topright)
plot!(tvec, cn, lab="DSM")
plot!(tvec, cn_static, lab="Static")
# display(cnplt) 


cyclecnplt = plot(xaxis="Angle of Attack (deg)", yaxis=L"C_l", leg=:topleft)
plot!((alphavec.*(180/pi)), cn, lab="DSM's Øye", linewidth=2)
plot!(alphavec.*(180/pi), cn_static, lab="Static", linewidth=0.5)
#display(cyclecnplt) 

#=
I think that I have the models matching, the difference is that I'm using Cn in this model and Cl in the other. 
=#

Faber = readdlm("Faber_Oye_0015_Full.csv" , ',')

x = Faber[:,1]
y = Faber[:,2]

Experimental = readdlm("Experimental_0015_Full.csv" , ',')

x1 = Experimental[:,1]
y1 = Experimental[:,2]

plot!(x1,y1, linestyle=:dash, linewidth = 2, label = "Experimental")

full_plot = scatter!(x,y, label = "Faber's Øye")

display(full_plot)

#savefig("full_plot_A_8")
"""
dsm.separationpoint.(Ref(af), polar_0015[52:end,1])
thing = plot(polar_0015[52:end,1].*180/pi , dsm.separationpoint.(Ref(af), polar_0015[52:end,1]))

display(thing)
"""

nothing