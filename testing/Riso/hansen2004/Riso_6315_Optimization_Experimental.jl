path = dirname(@__FILE__)
cd(path)
dsm = DynamicStallModels



file2 = "../../../polars/Hansen_6315_Top_Half_Blue.csv"
Hansen_Data = readdlm(file2, ',')
#Hansen_Data[:,1] = Hansen_Data[:,1].*pi/180


file = "../../../polars/Hansen_6315.csv"
polar = readdlm(file , ',')


c = 0.55
M = 0.1
v = 343
Vrel = M*v

tspan = (0.0, 2.0) 

answer = zeros(25, 2)
Error = zeros(1, length(Hansen_Data[:,1]))

dsmodel = Riso(Functional(), [0.294, 0.331] , [0.1, 0.3], [1.5, 6.0]) 

airfoils = Array{Airfoil, 1}(undef, 1)


function Uvector(t)
    return 343*0.1
end


parameters = [Uvector, 0.0, alphavec, alphavecdot]

x_initial = [0.0, 0.0, 0.0, 0.0]


function objective(g, x)

    x1, x2, x3, x4, x5, x6, x7 = x


    
    dsmodel.A[1] = x1
    dsmodel.A[2] = x2
    dsmodel.b[1] = x3
    dsmodel.b[2] = x4
    dsmodel.T[1] = x5
    dsmodel.T[2] = x6



    af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.RSP())


    airfoils[1] = af

    function alphavec(t)
        c = x7
        M = 0.1
        a = 343.0
        shift = 12.0
        amp = 4.0
        k = 0.1
    
        v = M*a
        omega = k*2*v/c
    
        alf = shift + amp*sin(omega*t)
        return alf*(pi/180)
    end
    
    function alphavecdot(t)
        c = x7
        M = 0.1
        a = 343.0
        shift = 12.0
        amp = 4.0
        k = 0.1
    
        v = M*a
        omega = k*2*v/c
    
        alf = amp*omega*cos(omega*t)
        return alf*(pi/180)
    end


    prob = ODEProblem(airfoils, x_initial, tspan, parameters)
    sol = DifferentialEquations.solve(prob, reltol=1e-8)

    #answer = parsesolution(dsmodel, airfoils, sol, parameters)

    dcldalpha = af.dcldalpha
    alpha0 = af.alpha0
    A1 = x1
    A2 = x2
    n = 1

    for i in 0.89:0.01:1.13
        ae = alphavec(i)*(1-A1-A2) + sol(i)[1] + sol(i)[2]
        Cl_fs = dsm.hansen_fully_sep(af, ae)
        Tu = c/(Uvector(i)*2)

        Cl_Dyn = dcldalpha*(ae - alpha0)*sol(i)[4] + Cl_fs*(1 - sol(i)[4]) + pi*Tu*alphavecdot(i)
        answer[n, 1] = alphavec(i)
        answer[n, 2] = Cl_Dyn

        
        n = n + 1
    end


    cl_fit = Akima(answer[:,1], answer[:,2])



    for i in 1:length(Hansen_Data[:,1])
        Error[1,i] = abs(Hansen_Data[i,2] - cl_fit(Hansen_Data[i,1]*pi/180))
    end

    g[1] = x3 - x4
    

    return sum(Error)

end

