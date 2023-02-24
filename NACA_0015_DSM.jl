using DelimitedFiles

cd("/Users/westonpace/Desktop")

df = readdlm("NACA_0015.csv" , ',')

matrix = []


for i in 1:(size(df,1))
    push!(matrix, [(pi/180)*df[i,1], df[i,2], df[i,3], 0.0])
end

polar_0015 = reduce(vcat,transpose.(matrix))

println(polar_0015)

println(length(polar_0015))