using DelimitedFiles

cd("/Users/westonpace/Desktop")

df = readdlm("Larsen_Data.csv" , ',', skipstart=1)

matrix = []

for i in 1:size(df,1)
    push!(matrix, [(pi/180)*df[i,1], df[i,2], 0.0, 0.0])
end

println(matrix)
