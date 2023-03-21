using DelimitedFiles, Plots, FLOWMath

cd("/Users/westonpace/Desktop/Oye_Error")

computed = readdlm("Faber_0015_Upstroke.csv" , ',')

x = computed[:,1]
y = computed[:,2]

Upstroke = Akima(x,y)

Faber = readdlm("Faber_Oye_0015.csv" , ',')

x1 = Faber[:,1]
y1 = Faber[:,2]

SE_Computed = []

for i in 1:length(x1)
    push!(SE_Computed, abs(y1[i]-Upstroke(x1[i])))
end

#abs_residual = sum(SE_Computed)
 
#abs_residual in the upstroke region after linear = 1.6584

Experimental = readdlm("Experimental_0015_Upstroke.csv" , ',')

x2 = Experimental[:,1]
y2 = Experimental[:,2]

Experimental_Fit = Akima(x2,y2)

SE_Experimental = []

for i in 1:length(x1)
    push!(SE_Experimental, abs(Experimental_Fit(x1[i])-Upstroke(x1[i])))
end

#println(sum(SE_Experimental))

#SE_Experimental (absolute residual between computed and experimental) for upstroke is = 1.8714

SE_Experimental_Faber = []

for i in 1:length(x1)
    push!(SE_Experimental_Faber, abs(Experimental_Fit(x1[i])-y1[i]))
end

println(sum(SE_Experimental_Faber))

#SE_Experimental_Faber (absolute residual between Faber's and Expeerimental) for upstroke is = 3.5112

"""
I have only done the residual for value in the post linear upstroke region; I will do the down stroke region later.
"""



