using Plots, Statistics, DelimitedFiles, Roots

include("../../Riso.jl")

k = 0.1 #given
u = 1.0 #Assumed
# u = 10.9236 #Optimized 2
c = 1.0 #Assumed
# c = 1.1474 #Optimized
# c = 1.07942 #Optimized 2
#opt 3
# u =  1.029259
# c = 0.54307

omega = k*2*u/c #implied

function U(t)
    return u
end

function Udot(t)
    return 0.0
end

function V(t)
    return 0.0
end

function alpha(t)
    # k = 0.1 #given
    # u = 50 #Assumed
    # c = 1.0 #Assumed

    # omega = k*2*u/c #implied

    amp = 4
    shift = 12

    alpha = amp*sin(omega*t) + shift
    return alpha*(pi/180)
end

function alphadot(t)
    # k = 0.1 #given
    # u = 50 #Assumed
    # c = 1.0 #Assumed

    # omega = k*2*u/c #implied

    amp = 4
    shift = 12

    alphadot = amp*omega*cos(omega*t)

    return alphadot*(pi/180)
end

polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/polar.csv", ',')

polar[:,1] = polar[:,1].*(pi/180)
liftfit = Akima(polar[:,1], polar[:,2])

dcldalpha = 2*pi*1.05
alpha0 = find_zero(liftfit, 0.0)

yintercept = -alpha0*dcldalpha

linearlift(alpha) = dcldalpha*alpha + yintercept
linearcl = linearlift.(polar[:,1])
staticlift = liftfit.(polar[:,1])

staticplt = plot(legend=:topleft, title="Static Lift Plot") #, xlim=(-0.07,0.25))
scatter!(polar[:,1], polar[:,2], lab="Static")
plot!(polar[:,1], linearcl, lab="Linear")
plot!(polar[:,1], staticlift, lab="liftfit")
# display(staticplt)


# A = [0.294, 0.331] #From the Hansen 2004 paper, for Riso airfoil
# b = [0.0664, 0.3266] # "" ""
# # Tp = 1/0.4125 #From the Larsen 2007 paper for the Vertol airfoil
# # Tf = 1/0.0875 # "" "" 
# # A = [ 0.3, 0.7] #Leishman 1990
# # b = [0.14, 0.53] # "" "" 
# Tp = 3.0 # "" "" 
# Tf = 6.0 # "" "" 

#Optimized coefficients #It caught the incoming condition. Maybe I need to make sure that it has converged before it goes crazy. 
# A = [0.9648, 0.34807]
# b = [0.07895, 0.2035]
# Tp = 2.988
# Tf = 5.98s
#Second optimized coefficients # A similar thing happened. 
# A = [0.256919, 0.2998237]
# b = [0.0835, 0.097929]
# Tp = 3.00979
# Tf = 5.9185
#Opt 3
# A = [0.28433, 0.7157]
# b = [0.01, 0.03884]
# Tp = 2.988
# Tf = 5.98
#Opt 4
# A = [0.3354261839159328, 0.3356162905032279]
# b = [0.010121096168090178, 0.2]
# Tp = 0.5
# Tf = 0.7932678364562898
#Opt 5
A = [0.1, 0.1]
b = [0.018617983848921067, 0.2]
Tp = 1.0
Tf = 1.0


x0 = zeros(4)
p = [c, A, b, Tp, Tf, dcldalpha, liftfit, U, Udot, V, alpha, alphadot, alpha0]
tspan = (0.0, 80.0)

prob = ODEProblem(states!, x0, tspan, p)
sol = solve(prob, dtmax=0.1)

Cld1, Cdd1, Cmd1, u1, t1 = parsesolution(sol, p, polar)
alphavec = alpha.(sol.t)
alfavec = alphavec.*(180/pi)



# statesplt = plot(sol,linewidth=2,xaxis="t",label=["x1" "x2" "x3" "x4"],layout=(4,1))
# display(statesplt)


colors = theme_palette(:auto).colors.colors


explift = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/circl12pm4.csv", ',')
tidx = 400

clplt = plot(legend=:topleft, yaxis="Cl", xaxis="Alpha (deg)", xlims=(5.0, 20.0), ylims=(0.75, 1.75))
plot!(polar[:,1].*(180/pi), polar[:,2], lab="Static Lift", linestyle=:dot, linewidth=4, color=colors[3])
scatter!(explift[:,1], explift[:,2], lab="Hansen 2004", markersize=4.5, markeralpha=0.9, color=colors[1])
plot!(alfavec[tidx:end], Cld1[tidx:end], lab="Present", linewidth=3, color=colors[2])
display(clplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/riso/nonlinearsolve/stalllift_opt.png")

expdrag = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/drag_12pm4.csv", ',')

cdplt = plot(legend=:topleft, yaxis="Cd", xaxis="Alpha (deg)", xlims=(5.0, 20.0), ylims=(-0.01, 0.16))
plot!(polar[:,1].*(180/pi), polar[:,3], lab="Static Drag", color=colors[3], linestyle=:dot, linewidth=4)
scatter!(expdrag[:,1], expdrag[:,2], lab="Hansen 2004", color=colors[1], markersize=4.5, markeralpha=0.9)
plot!(alfavec[tidx:end], Cdd1[tidx:end], lab="Present", linewidth=3, color=colors[2])
display(cdplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/riso/nonlinearsolve/stalldrag_opt.png")

expmoment = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/moment12pm4.csv", ',')

cmplt = plot(legend=:topleft, title="Cyclic Moment", yaxis="Cm", xaxis="Alpha (deg)", xlims=(-5.0, 30))
scatter!(expmoment[:,1], expmoment[:,2], lab="Paper Values")
plot!(alfavec, Cmd1, lab="My Values")
plot!(polar[:,1].*(180/pi), polar[:,4], lab="Static Moment")
# display(cmplt)
# savefig("/Users/adamcardoza/Box/research/FLOW/bladeopt/figures/dynamicstall/riso/nonlinearsolve/stallmoment_opt.png")



nothing