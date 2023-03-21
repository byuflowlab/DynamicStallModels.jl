using DelimitedFiles
using CCBlade
using Plots
using FLOWMath

cd("/Users/westonpace/Desktop")

df = readdlm("NACA_0015_Best_Fit.csv" , ',')

af =readdlm("NACA_0015_Best_Fit_Drag.csv", ',')


alpha_0 = df[:,1]
cl_0 = df[:,2]
alpha_1 = af[:,1]
cd_0 = af[:,2]
cr75 = 1000


cl_fit = Akima(alpha_0, cl_0)
cd_fit = Akima(alpha_1, cd_0)

a_vec = 0.0:0.25:20

cl_fit_values = cl_fit.(a_vec)

cd_fit_values = cd_fit.(a_vec)

aoa, lift, drag = viterna(a_vec.*pi/180, cl_fit_values, cd_fit_values, cr75)

matrix = []

for i in 1:length(aoa)
    push!(matrix, [aoa[i], lift[i], drag[i], 0.0])
end

polar_0015 = reduce(vcat,transpose.(matrix))

plot(polar_0015[:,1], polar_0015[:,3])