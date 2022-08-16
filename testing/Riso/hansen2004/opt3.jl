using Plots, Statistics, DelimitedFiles, Roots, Optim, Snopt, FiniteDiff

include("../Riso.jl")

k = 0.1 #given
u = 1.0
c = 1.0

omega = k*2*u/c #implied

polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/polar.csv", ',')
    
polar[:,1] = polar[:,1].*(pi/180)
liftfit = Akima(polar[:,1], polar[:,2])

dcldalpha = 2*pi*1.05
alpha0 = find_zero(liftfit, 0.0)

yintercept = -alpha0*dcldalpha

linearlift(alpha) = dcldalpha*alpha + yintercept
linearcl = linearlift.(polar[:,1])
staticlift = liftfit.(polar[:,1])

explift = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/circl12pm4.csv", ',')

expdrag = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/drag_12pm4.csv", ',')

expmoment = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/moment12pm4.csv", ',')

farval, fariter = findmax(explift[:,1]) 
alpha1 = farval
t1 = (asin((alpha1-12)/4)+ 4*pi)/omega

closeval, closeiter = findmin(explift[:,1])
alpha2 = closeval
t2 = (asin((alpha2-12)/4)+ 4*pi)/omega

maxval, maxiter = findmax(explift[:,2])
alpha3 = explift[maxiter,1]
t3 = (asin((alpha3-12)/4)+ 4*pi)/omega

minval, miniter = findmin(explift[:,2])
alpha4 = explift[miniter,1]
t4 = (asin((alpha4-12)/4)+ 4*pi)/omega

function getobj(x)
    
    A1, A2, b1, b2, Tp, Tf = x

    A = [A1, A2]
    b = [b1, b2]

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

    Cld1, Cdd1, Cmd1, u1, t1vec = parsesolution(sol, p, polar)
    # println(typeof(alphafun))
    alphavec = alphafun.(t1vec)
    alfavec = alphavec.*(180/pi)

    index1 = findall(x->isapprox(x,t1,atol=1e-1),sol.t)[1]
    dL1 = abs(explift[fariter,2]-Cld1[index1])
    dD1 = abs(expdrag[fariter,2]-Cdd1[index1])
    dM1 = abs(expmoment[fariter,2]-Cmd1[index1])

    index2 = findall(x->isapprox(x,t2,atol=1e-1),sol.t)[1]
    dL2 = abs(explift[closeiter, 2]-Cld1[index2])
    dD2 = abs(expdrag[closeiter,2]-Cdd1[index2])
    dM2 = abs(expmoment[closeiter,2]-Cmd1[index2])

    index3 = findall(x->isapprox(x,t3,atol=1e-1),sol.t)[1]
    dL3 = abs(explift[maxiter, 2]-Cld1[index3])
    dD3 = abs(expdrag[maxiter,2]-Cdd1[index3])
    dM3 = abs(expmoment[maxiter,2]-Cmd1[index3])

    index4 = findall(x->isapprox(x,t4,atol=1e-1),sol.t)[1]
    dL4 = abs(explift[miniter, 2]-Cld1[index4])
    dD4 = abs(expdrag[miniter,2]-Cdd1[index4])
    dM4 = abs(expmoment[miniter,2]-Cmd1[index4])

    davg = dL1 + dL2 + dD1 + dD2 + dM1 + dM2 + dL3 + dD3 + dM3 + dL4 + dD4 + dM4
    return davg
end

println("Preparing Optimization")
A = [0.294, 0.331] #From the Hansen 2004 paper, for Riso airfoil
b = [0.0664, 0.3266] # "" ""
Tp = 3.0 # "" "" 
Tf = 6.0 # "" "" 

X0 = [A[1], A[2], b[1], b[2], Tp, Tf]
XL = [0.1, 0.1, 0.01, 0.2, 1.0, 1.0]
XU = [1.0, 1.0, 0.5, 0.5, 30.0, 30.0]

function confunc(x)
    c = zeros(2)
    c[1] = x[1] + x[2] - 1.0 #x[1] and x[2] must be less than 1.0. 
    c[2] = x[1] - x[2] #x[1] must be less than x[2]
    return c 
end

function wrapper!(y,x)
    y[1] = getobj(x)
    y[2:end] = confunc(x)
end

function objcon(x)
    f = getobj(x)
    g = confunc(x)
    FiniteDiff.finite_difference_jacobian!(J, wrapper!, x; relstep=0.1)

    dfdx = J[1, :]
    dgdx = J[2:end, :]
    
    fail = false
    return f, g, dfdx, dgdx, fail
end

J = zeros(3, length(X0))

##SNOPT options
options = Dict{String, Any}()
options["Derivative option"] = 1
options["Verify level"] = 1
options["Print file"] = "printfile.txt"
options["Summary file"] = "sumfile.txt"
### Run the optimization
println("Optimizing...")

xopt, fopt, info = snopt(objcon, X0, XL, XU, options)


nothing