using DelimitedFiles, Plots, Statistics, FLOWMath, DifferentialEquations, Roots

#=
Validate the Riso model's ability to simulate experimental data from Larsen's 2007 paper. 


Adam Cardoza
7/27/22
=#

include("../riso.jl")
include("../indicialriso.jl")

# expdata = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/experimentaldata/Larsen2007/riso1.csv", ',') #Todo: Move this to a relative reference (inside the repo). 

expdata = readdlm("../../../experimentaldata/Larsen2007/Riso/stall/riso1.csv", ',')

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

function alpha(t)
    k = 0.062
    v = 60
    c = 1.5

    amp = 4.91
    shift = 7.33

    omega = 2*k*v/c

    alfa =  amp*cos(omega*t) + shift
    return alfa*pi/180
end

function alphadot(t)
    k = 0.062
    v = 60
    c = 1.5

    amp = 4.91
    shift = 7.33

    omega = 2*k*v/c

    alfadot =  -amp*sin(omega*t)/omega
    return alfadot*pi/180
end


#Enviroment
v = 60

#Geometry
c = 1.5
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/learning/exampledata/xf-v23010-il-200000-n5.csv", ','; skipstart=12)
# liftfit = Akima(polar[:,1].*(pi/180), polar[:,2])
# polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/polars/extendedVertol 23010-1.58.dat", ',')
polar = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/data/polars/extendedVertol 23010-1.58.dat", ',') #Todo: Polar 1 is a actually a little off from the experimental data in Larsen's paper. 
# polar[:,1] = polar[:,1].*(pi/180)

polar = hcat(polar, zeros(length(polar[:,1])), zeros(length(polar[:,1])))

liftfit = Akima(polar[:,1], polar[:,2])
dcldalpha = 2*pi*1.05
alpha0 = 0.014 # -0.019

alphas = [find_seperation_alpha(liftfit, dcldalpha, alpha0)...]
afm = alphas[2]
afp = alphas[1]

### Airfoil Constants
A = [0.165, 0.335] #Todo: I haven't checked these out. -> They seem to be doing well. 
b = [0.0455, 0.3] 
Tp = 1/0.4125
Tf = 1/0.0875


### Initialize 
x0 = zeros(4)
x0[4] = 0.98
p = [U, Udot, alpha, alphadot, c, A[1], A[2], b[1], b[2], dcldalpha, alpha0, Tp, Tf, liftfit, afm, afp]
tspan = (0.0, 20.0)

prob = ODEProblem(riso_ode!,x0,tspan,p)
sol = solve(prob)


### Solve using the indicial method #Todo: 
states = indicialsolve(sol.t, U, Udot, alpha, alphadot, c, dcldalpha, alpha0, afm, afp, A, b, Tp, Tf, liftfit, x0)

########### Parse the solution ############
Cld, Cdd, Cmd, u, t = parsesolution(sol, p, polar)
Cldi, Cddi, Cmdi = parseindicialsolution(states, p, polar)



alphavec = alpha.(sol.t)
Cls = liftfit.(alphavec)
alfavec = alphavec.*(180/pi)
staticcl = liftfit.(alphavec)


######### Plotting #############

### Plot the states
statesplt = plot(sol,linewidth=2,xaxis="t",label=["x1" "x2" "x3" "x4"],layout=(4,1))
plot!(sol.t, states)
# display(statesplt)

### Plot the lift as a function of time. 
clplt = plot(sol.t, Cld, lab="dynamic", leg=:bottomright)
plot!(sol.t, Cls, lab="static") 
xlabel!("time (s)")
ylabel!("Coefficient of Lift")
# display(clplt)



### Plot the lift as a function of angle of attack. 
cycleplt = plot(legend=:topleft)
plot!(alfavec, Cld, lab="Riso")
scatter!(expdata[:,1], expdata[:,2], lab="Experimental")
plot!(alfavec, staticcl, lab="Static")
plot!(alfavec, Cldi, lab="Indicial Riso")
display(cycleplt)


nothing