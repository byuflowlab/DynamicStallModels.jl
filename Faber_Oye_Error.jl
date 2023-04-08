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

#println(sum(SE_Experimental_Faber))

#SE_Experimental_Faber (absolute residual between Faber's and Experimental) for upstroke is = 3.5112

"""
I have only done the residual for value in the post linear upstroke region; I will do the down stroke region later.
"""

computed_Downstroke = readdlm("Faber_0015_Downstroke.csv" , ',')

x_d = computed_Downstroke[:,1]
y_d = computed_Downstroke[:,2]

Downstroke = Akima(x_d,y_d)

Faber_Downstroke = readdlm("Faber_Oye_0015_Downstroke.csv" , ',')

x1_D = Faber_Downstroke[:,1]
y1_D = Faber_Downstroke[:,2]

SE_Computed_Downstroke = []

for i in 1:length(x1_D)
    push!(SE_Computed_Downstroke, abs(y1_D[i]-Downstroke(x1_D[i])))
end

#println(sum(SE_Computed_Downstroke))
#residual between the computed and faber's oye in the downstroke region = 0.88799

Experimental_Downstroke = readdlm("Experimental_0015_Downstroke.csv" , ',')

x2_D = Experimental_Downstroke[:,1]
y2_D = Experimental_Downstroke[:,2]

Experimental_Fit_Downstroke = Akima(x2_D,y2_D)

SE_Experimental_Downstroke = []

for i in 1:length(x1_D)
    push!(SE_Experimental_Downstroke, abs(Experimental_Fit_Downstroke(x1_D[i])-Downstroke(x1_D[i])))
end

#println(sum(SE_Experimental_Downstroke))
#the residual between the experimental and computed in the downstroke region is = 3.5606

SE_Experimental_Faber_Downstroke = []

for i in 1:length(x1_D)
    push!(SE_Experimental_Faber_Downstroke, abs(Experimental_Fit_Downstroke(x1_D[i])-y1_D[i]))
end

#println(sum(SE_Experimental_Faber_Downstroke))
#residual between Faber's oye and thee experimental in the downstroke = 3.5998