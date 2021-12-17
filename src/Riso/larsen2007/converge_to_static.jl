using DelimitedFiles, Plots, Statistics

include("../../Riso.jl")

expdata = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Larsen2007/riso1.csv", ',')

middle = (maximum(expdata[:,1])+minimum(expdata[:,1]))/2
delalpha = maximum(expdata[:,1])-middle


function U(t)
    return 60
end

function Udot(t)
    return 0
end

function V(t)
    return 0.0
end

function alphadot(t)

    return 0
end


#Enviroment
v = 60

#Geometry
c = 1.5
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/learning/exampledata/xf-v23010-il-200000-n5.csv", ','; skipstart=12)
# liftfit = Akima(polar[:,1].*(pi/180), polar[:,2])
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/extendedVertol 23010-1.58.dat", ',')
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Larsen2007/vertol_lowaoa_static.csv", ',') 
polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Larsen2007/fullpolar.csv", ',') 
# polar[:,1] = polar[:,1].*(pi/180)

polar = hcat(polar, zeros(length(polar[:,1])), zeros(length(polar[:,1])))

# liftfit = Akima(polar[:,1], polar[:,2])
liftfit = Akima(polar[:,1], polar[:,2])
dcldalpha = 2*pi*1.05
alpha0 = 0.014 # -0.019

#Constants
A = [0.165, 0.335] 
b = [0.0455, 0.3] 
Tp = 1/0.4125
Tf = 1/0.0875

#Initialize 

tspan = (0.0, 300.0)

clvec = zeros(length(polar[:,1]))

for i = 1:length(clvec)
    
    alpha(t) = polar[i,1]

    x0 = zeros(4)
    p = [c, A, b, Tp, Tf, dcldalpha, liftfit, U, Udot, V, alpha, alphadot, alpha0]

    prob = ODEProblem(states!,x0,tspan,p)
    sol = solve(prob)

    Cld, Cdd, Cmd, u, t = parsesolution(sol, p, polar)
    if i==length(clvec)
        plt = plot(t, Cld)
        display(plt)
    end
    clvec[i] = Cld[end]
end

alphavec = polar[:,1].*(180/pi)

staticplt = plot(legend=:topleft)
plot!(alphavec, clvec, lab="Riso")
plot!(alphavec, polar[:,2], lab="Static")
display(staticplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/riso/larsen/convergedstatic_lift.png")


nothing