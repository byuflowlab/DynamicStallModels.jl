using Plots, Statistics, DelimitedFiles, Roots, DynamicStallModels, DifferentialEquations

DE = DifferentialEquations
ds = DynamicStallModels


k = 0.2 #given
u = 1.0 #Assumed
c = 1.0 #Assumed

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

    amp = 2
    shift = 5

    alpha = amp*sin(omega*t) + shift
    return alpha*(pi/180)
end

function alphadot(t)
    # k = 0.1 #given
    # u = 50 #Assumed
    # c = 1.0 #Assumed

    # omega = k*2*u/c #implied

    amp = 2
    shift = 5

    alphadot = amp*omega*cos(omega*t)

    return alphadot*(pi/180)
end

linearalpha = -pi:0.01:pi
polar = hcat(linearalpha, 2*pi.*linearalpha, zeros(length(linearalpha)), zeros(length(linearalpha)))

# polar = readdlm("./data/Hansen2004/figure8_flatplate/flatplatepolar.csv", ',')
# polar[:,1] .*= pi/180

A = [0.165, 0.335] #From the Hansen 2004 paper, for flat plate
b = [0.0455, 0.3000] # "" ""
# A = [ 0.3, 0.7] #Leishman 1990
# b = [0.14, 0.53] # "" "" 
Tp = 3.0 # "" "" 
Tf = 6.0 # "" "" 
T = [Tp, Tf]


m, n = size(polar)
newpolar = hcat(polar, zeros(m), zeros(m))
afs = Array{Airfoil, 1}(undef, 1)
afs[1] = airfoil(newpolar; A, b, T)

x0 = zeros(4)
x0[4] = 1.0


tspan = (0.0, 80.0)
dt = 0.01
t0 = tspan[1]

dsmodel = riso(afs)

p = [U(t0), Udot(t0), alpha(t0), alphadot(t0), c]


tvec = tspan[1]:dt:tspan[2]
nt = length(tvec)

prob = ODEProblem{false}(dsmodel, x0, tspan, p)
integrator = DE.init(prob, Tsit5(); verbose=false)

xds = zeros(nt, 4)
xds[1,:] = x0

cl_d = zeros(nt)
cd_d = zeros(nt)

cl_i, cd_i = ds.parsestates(dsmodel, xds[1,:], integrator.p)
cl_d[1] = cl_i[1]
cd_d[1] = cd_i[1]

for i = 2:nt
    t = tvec[i]
    dt_local = tvec[i]-tvec[i-1]
    integrator.p[1] = U(t)
    integrator.p[2] = Udot(t)
    integrator.p[3] = alpha(t)
    integrator.p[4] = alphadot(t)

    DifferentialEquations.step!(integrator, dt_local, true)
    xds[i,:] = integrator.u

    cl_local, cd_local = ds.parsestates(dsmodel, xds[i,:], integrator.p)

    cl_d[i] = cl_local[1]
    cd_d[i] = cd_local[1]
end



statesplt = plot(xaxis="Time (s)", yaxis="States")
plot!(tvec, xds)
# display(statesplt)


alphavec = alpha.(tvec)

expdata = readdlm("../../../data/Hansen2004/figure8_flatplate/indicial.csv", ',')
staticdata = dsmodel.airfoils[1].cl.(alphavec)


clplt = plot(legend=:topleft, title="Cyclic Alpha", yaxis="Cl", xaxis="Alpha (deg)")
scatter!(expdata[:,1], expdata[:,2], lab="Paper Values")
plot!(alphavec.*(180/pi), cl_d, lab="My Values")
plot!(alphavec.*(180/pi), staticdata, lab="Static")
display("image/png", clplt)

#=
#Todo. I need better static data for a good match. 
    -> Much better when I just use inviscid flow theory (m=2*pi).

=#


nothing