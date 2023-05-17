#=
Todo: What is this doing here? I feel like this file belongs in the testing section. 
=#

using DelimitedFiles, Plots, Aplots, Revise



bll = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/outputs/leishman1989/naca0012_8.1_10.3s_Cdn_BLL.csv", ',')
leishman = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/outputs/leishman1989/naca0012_8.1_10.3s_Cdn_BLss.csv", ',')
oye = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/outputs/leishman1989/naca0012_8.1_10.3s_Cdn_Oye.csv", ',')
larsen = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/outputs/leishman1989/naca0012_8.1_10.3s_Cdn_Larsen.csv", ',')
exppolar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989fig8Cn.csv", ',')

Cnplt = plot(leg=:topleft)
scatter!(exppolar[:,1], exppolar[:,2], lab="Experimental", markercolor=:red3)
aplot!(bll[100:end,1], bll[100:end,2], legend=:topleft, skip=true, skipnum=15, label="Leishman-Larsen", color=:dodgerblue2)
aplot!(leishman[100:end,1], leishman[100:end,2],legend=:topleft, label="Leishman", markershape=:utriangle, color=:green, skip=true, skipnum=15)
aplot!(oye[:,1], oye[:,2], legend=:topleft, label="Oye", markershape=:cross, color=:purple, marker=true)
aplot!(larsen[100:end,1], larsen[100:end,2], legend=:topleft, label="Larsen", skip=true, skipnum=15, markershape=:diamond, color=:tan1)
ylabel!("Normal Coefficient")
xlabel!("Angle of Attack (degrees)")
# display(Cnplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/comparison/exp_BLL_BLss_Oye_Larsen.png")

logocolors = Colors.JULIA_LOGO_COLORS

leishplt = plot(leg=false)
scatter!(exppolar[:,1], exppolar[:,2], lab="Experimental", markercolor=:orangered)
aplot!(leishman[100:end,1], leishman[100:end,2],legend=false, label="Leishman", markershape=:utriangle, color=:green, skip=true, skipnum=15)
ylims!((0.2,1.9))

BLLplt = plot(leg=false)
scatter!(exppolar[:,1], exppolar[:,2], lab="Experimental", markercolor=:orangered)
aplot!(bll[100:end,1], bll[100:end,2], legend=false, skip=true, skipnum=15, label="Leishman-Larsen", color=:dodgerblue2)
ylims!((0.2,1.9))

Oyeplt = plot(leg=false)
scatter!(exppolar[:,1], exppolar[:,2], lab="Experimental", markercolor=:orangered)
aplot!(oye[:,1], oye[:,2], legend=false, label="Oye", markershape=:cross, color=:purple, marker=true)
ylims!((0.2,1.9))

Larsenplt = plot(leg=false)
scatter!(exppolar[:,1], exppolar[:,2], lab="Experimental", markercolor=:orangered)
aplot!(larsen[100:end,1], larsen[100:end,2], legend=false, label="Larsen", skip=true, skipnum=15, markershape=:diamond, color=:tan1)
ylims!((0.2,1.9))

compplt = plot(leishplt, BLLplt, Larsenplt, Oyeplt, layout=(2,2), title=["Beddoes-Leishman" "Leishman-Larsen" "Larsen" "Oye"], size=(1200,800), xaxis="Angle of Attack (degrees)", yaxis="Coefficient of Lift")
display(compplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/comparison/exp_BLL_BLss_Oye_Larsen_sidebyside.png")



nothing

