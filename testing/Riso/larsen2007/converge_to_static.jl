using DelimitedFiles, Plots, Statistics, FLOWMath, DifferentialEquations, Roots

cd("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/dynamicstallmodels/src/Riso/larsen2007")

include("../riso.jl")

expdata = readdlm("../../../experimentaldata/Larsen2007/Riso/stall/riso_stall_experimental_Larsen2007.csv", ',')
paperdata = readdlm("../../../experimentaldata/Larsen2007/Riso/stall/riso_stall_paper_Larsen2007.csv", ',')



function U(t)
    return 60
end

function Udot(t)
    return 0
end

function alphadot(t)

    return 0
end


#Enviroment
v = 60

#Geometry
c = 1.5
polar = readdlm("../../../polars/extendedVertol 23010-1.58.dat", ',')
liftfit = Akima(polar[:,1], polar[:,2])
dragfit = Akima(polar[:,1], polar[:,3])
dcldalpha = 2*pi*1.05
alpha0 = -0.019
Cd0 = dragfit(alpha0)

Cli(alpha) = dcldalpha*(alpha-alpha0)

alphas = [find_seperation_alpha(liftfit, dcldalpha, alpha0)...]
afm = alphas[2]
afp = alphas[1]

#Constants
A = [0.165, 0.335] 
b = [0.0455, 0.3] 
Tp = 1/0.4125
Tf = 1/0.0875

#Initialize 

tspan = (0.0, 300.0)

aoavec = -50:0.5:50
alphavec = aoavec.*(pi/180)

nn = length(aoavec)
clvec = zeros(nn)

for i = 1:nn
    
    alpha(t) = alphavec[i]

    x0 = zeros(4)
    x0[4] = 0.98
    p = [U, Udot, alpha, alphadot, c, A[1], A[2], b[1], b[2], dcldalpha, alpha0, Tp, Tf, liftfit, afm, afp]
    

    prob = ODEProblem(riso_ode!,x0,tspan,p)
    sol = solve(prob)


    Cld, cdvec, t = parsesolution(sol, p, dragfit, Cd0)

    if i==nn
        # @show alpha(0)
        # plt = plot(t, Cld)
        # display(plt)
        # @show sol[end]
        plt2 = plot(sol)
        display(plt2)
    end
    clvec[i] = Cld[end]
end

clivec = Cli.(alphavec)


### Test the final point to find what's going on in the solution. 
# alpha2(t) = alphavec[end]
# x02 = zeros(4)

# #Todo: The fourth state should descend to zero, because we are well in the stall region. 
# # x02[4] = 0.98
# p2 = [U, Udot, alpha2, alphadot, c, A[1], A[2], b[1], b[2], dcldalpha, alpha0, Tp, Tf, liftfit, afm, afp]

# prob2 = ODEProblem(riso_ode!,x02,tspan,p2)
# sol2 = solve(prob2)

# Cld2, cdvec2, t2 = parsesolution(sol2, p2, dragfit, Cd0)

staticplt = plot(legend=:topleft, xlim=(-50,50), ylim=(-1.0,1.5))
plot!(aoavec, clvec, lab="Riso")
plot!(polar[:,1].*(180/pi), polar[:,2], lab="Static")
plot!(aoavec, clivec, lab="Inviscid")
vline!(alphas.*(180/pi), lab="Seperation Angles")
display(staticplt)



nothing