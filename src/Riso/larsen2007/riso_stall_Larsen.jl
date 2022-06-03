using DifferentialEquations, Plots, DelimitedFiles, FLOWMath, Roots


include("../riso.jl")

#=
Use the attached region experimental data and results from Larsen's 2007 paper, to replicate results and validate my model. 
=#

expdata = readdlm("../../../experimentaldata/Larsen2007/Riso/stall/riso_stall_experimental_Larsen2007.csv", ',')
paperdata = readdlm("../../../experimentaldata/Larsen2007/Riso/stall/riso_stall_paper_Larsen2007.csv", ',')



U(t) = 60
Udot(t) = 0.0

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
V = 60

#Geometry
c = 1.5
polar = readdlm("../../../polars/extendedVertol 23010-1.58.dat", ',')
liftfit = Akima(polar[:,1], polar[:,2])
dragfit = Akima(polar[:,1], polar[:,3])
dcldalpha = 2*pi*1.05
alpha0 = 0.014 # -0.019
Cd0 = dragfit(alpha0)

alphas = [find_seperation_alpha(liftfit, dcldalpha, alpha0)...]
afm = alphas[2]
afp = alphas[1]

#Constants
A = [0.165, 0.335] 
b = [0.0455, 0.3]
Tp = 1/0.4125
Tf = 1/0.0875


#Initialize 
x0 = zeros(4)
x0[4] = 0.98
p = [U, Udot, alpha, alphadot, c, A[1], A[2], b[1], b[2], dcldalpha, alpha0, Tp, Tf, liftfit, afm, afp]
tspan = (0.0, 20.0)

prob = ODEProblem(riso_ode!,x0,tspan,p)
sol = solve(prob)


clvec, cdvec, t = parsesolution(sol, p, dragfit, Cd0)

clplt = plot(t, clvec, xaxis="Time (s)", yaxis="Coefficient of Lift")
display(clplt)

alphavec = alpha.(t)
aoavec = alphavec.*(180/pi)
clcycle = plot(xaxis="Angle of Attack (rads)", yaxis="Coefficient of Lift", legend=:topleft)
plot!(aoavec, clvec, lab="Current")
scatter!(expdata[:,1], expdata[:,2], lab="Experimental") 
scatter!(paperdata[:,1], paperdata[:,2], lab="Larsen 2007") 
display(clcycle)


#=
This matches pretty well. It could match a little better. I'm not sure if it is something with my coefficients, or the method. 
=#


nothing