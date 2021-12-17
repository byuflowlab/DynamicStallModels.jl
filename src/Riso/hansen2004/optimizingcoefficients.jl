using Plots, Statistics, DelimitedFiles, Roots, Optim

include("../Riso.jl")

k = 0.1 #given

polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/static.csv", ',')
    
plr = deepcopy(polar)
plr[:,1] = plr[:,1].*(pi/180)
liftfit = Akima(plr[:,1], plr[:,2])

dcldalpha = 2*pi*1.05
alpha0 = find_zero(liftfit, 0.0)

yintercept = -alpha0*dcldalpha

linearlift(alpha) = dcldalpha*alpha + yintercept
linearcl = linearlift.(plr[:,1])
staticlift = liftfit.(plr[:,1])

expdata = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/circl12pm4.csv", ',')

function getobj(x)
    
    u, c, A1, A2, b1, b2, Tp, Tf = x

    A = [A1, A2]
    b = [b1, b2]

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

    function alphafun(t)
        amp = 4
        shift = 12

        alpha = amp*sin(omega*t) + shift
        return alpha*(pi/180)
    end

    function alphadotfun(t)
        amp = 4
        shift = 12

        alphadot = amp*omega*cos(omega*t)

        return alphadot*(pi/180)
    end

    x0 = zeros(4)
    p = [c, A, b, Tp, Tf, dcldalpha, liftfit, U, Udot, V, alphafun, alphadotfun, alpha0]
    tspan = (0.0, 80.0)

    prob = ODEProblem(states!, x0, tspan, p)
    sol = solve(prob, dtmax=0.1)

    Cld1, u1, f1, t1vec = parsesolution(sol, p)
    # println(typeof(alphafun))
    alphavec = alphafun.(t1vec)
    alfavec = alphavec.*(180/pi)

    farval, fariter = findmax(expdata[:,1])
    alpha1 = farval
    t1 = (asin((alpha1-12)/4)+ 2*pi)/omega
    index1 = findall(x->isapprox(x,t1,atol=1e-1),sol.t)[1]
    d1 = abs(expdata[fariter,2]-Cld1[index1])

    closeval, closeiter = findmin(expdata[:,1])
    alpha2 = closeval
    t2 = (asin((alpha2-12)/4)+ 2*pi)/omega
    index2 = findall(x->isapprox(x,t2,atol=1e-1),sol.t)[1]
    d2 = abs(expdata[closeiter, 2]-Cld1[index2])

    davg = (d1+d2)/2
    return davg
end

u = 1.1 #Assumed
c = 1.0 #Assumed
A = [0.294, 0.331] #From the Hansen 2004 paper, for Riso airfoil
b = [0.0664, 0.3266] # "" ""
Tp = 3.0 # "" "" 
Tf = 6.0 # "" "" 

X0 = [u, c, A[1], A[2], b[1], b[2], Tp, Tf]
XL = [1.0, 0.05, 0.1, 0.1, 0.01, 0.01, 0.5, 0.5]
XU = [100.0, 5.0, 1.0, 1.0, 0.5, 0.5, 30.0, 30.0]

inneropt = GradientDescent()
result = optimize(getobj, XL, XU, X0, Fminbox(inneropt))




nothing