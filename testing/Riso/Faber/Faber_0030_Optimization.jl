path = dirname(@__FILE__)
cd(path)
dsm = DynamicStallModels


file = "../../../polars/Extended_NACA_0030.csv"
polar = readdlm(file, ',')


file2 = "../../../polars/Faber_0030_Riso_Top_Half.csv"
Faber_Result = readdlm(file2, ',')


c = 0.55
M = 0.11
v = 343


tspan = (0.0, 2.0)
dsmodel = Riso(Functional(), [0.294, 0.331] , [0.0664, 0.3266], [1.5, 6]) 


airfoils = Array{Airfoil, 1}(undef, 1)


function Uvector(t)
    return 0.11*343
end


function alphavec(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.64
    k = 0.051

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end


function alphavecdot(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.64
    k = 0.051

    v = M*a
    omega = k*2*v/c

    alf = amp*omega*cos(omega*t)
    return alf*(pi/180)
end


parameters = [Uvector, 0.0, alphavec, alphavecdot]
x_initial = [0.0, 0.0, 0.0, 0.0]

answer = zeros(45, 2)
Error = zeros(1, length(Faber_Result[:,1]))

function objective(g,x)

    dsmodel.A[1] = x[1]
    dsmodel.A[2] = x[2]
    dsmodel.b[1] = x[3]
    dsmodel.b[2] = x[4]
    dsmodel.T[1] = x[5]
    dsmodel.T[2] = x[6]

    af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.RSP())

    airfoils[1] = af

    prob = ODEProblem(airfoils, x_initial, tspan, parameters)
    sol = DifferentialEquations.solve(prob, reltol=1e-8)

    dcldalpha = af.dcldalpha
    alpha0 = af.alpha0
    A1 = x[1]
    A2 = x[2]
    n = 1

    for i in 0.68:0.01:1.12
        ae = alphavec(i)*(1-A1-A2) + sol(i)[1] + sol(i)[2]
        Cl_fs = dsm.hansen_fully_sep(af, ae)
        Tu = c/(Uvector(i)*2)

        Cl_Dyn = dcldalpha*(ae - alpha0)*sol(i)[4] + Cl_fs*(1 - sol(i)[4]) + pi*Tu*alphavecdot(i)
        answer[n, 1] = alphavec(i)
        answer[n, 2] = Cl_Dyn

        
        n = n + 1
    end

    cl_fit = Akima(answer[:,1], answer[:,2])

    for i in 1:length(Faber_Result[:,1])
        Error[1,i] = abs(Faber_Result[i,2] - cl_fit(Faber_Result[i,1]*pi/180))
    end

    g[1] = x[3] - x[4]
    

    return sum(Error)
end

