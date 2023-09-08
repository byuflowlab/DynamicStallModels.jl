using DynamicStallModels, DelimitedFiles, FLOWMath, DifferentialEquations, Test, Plots

dsm = DynamicStallModels


path = dirname(@__FILE__)
cd(path)



@testset "Riso Separation Point" begin
    file = "../data/Hansen_Figure3.csv"

    polar = readdlm(file , ',')
    
    c = 0.55
    M = 0.11
    v = 343
    Vrel = M*v
    
    dsmodel = Riso(Functional(), [0.294, 0.331] , [0.0664, 0.3266], [1.5, 6]) 
    af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.RSP())
    
    airfoils = Array{Airfoil, 1}(undef, 1)
    airfoils[1] = af
    
    tspan = (0.0, 0.224)
    
    function Uvector(t)
        return 0.11*343
    end
    
    function alphavec(t)
        c = 0.55
        M = 0.11
        a = 343.0
        shift = -3.0
        amp = 46.0
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
        shift = -3.0
        amp = 46.0
        k = 0.051
    
        v = M*a
        omega = k*2*v/c
    
        alf = amp*omega*cos(omega*t)
        return alf*(pi/180)
    end
    
    function separationpoint(airfoil::Airfoil, alpha)
        if isapprox(alpha, -pi, atol=1e-4)
            println("Riso sep function called. ")
        end
    
    
        afm, _ = airfoil.alphasep
        clfit = airfoil.cl
        dcldalpha = airfoil.dcldalpha
        alpha0 = airfoil.alpha0
    
        f(x) = clfit(x) - dcldalpha*(x - alpha0)/4
    
        afp , _ = brent(f , 0.17 , 0.8726)
    
    
        if !(afm < alpha < afp) #Check if alpha is in the bounds. 
            # println("f was set to zero. ")
            return typeof(alpha)(0)
        end
        
        cl_static = clfit(alpha)
        cl_linear = dcldalpha*(alpha-alpha0)
        f = (2*sqrt(abs(cl_static/cl_linear))-1)^2
    
        
        if f>1 #Question: What if I don't return this? I might get Inf.... or possibly NaN... but I will less likely get 1.0... which is my problem child in the seperated coefficient of lift function. -> I fixed the fully seperated coefficient of lift function... I just plugged this function inside the other and simplified. 
            return typeof(alpha)(1.0)
        elseif isnan(f)
            # println("f return NaN")
            return typeof(alpha)(1.0)
        end
        
    
        return f
    end
    
    sep_vec = zeros(1, 70)
    angles = zeros(1,70)
    n = [1]
    
    for i in 0:0.01:0.69
        sep_vec[1, n[1]] = separationpoint(af, i)
        angles[1, n[1]] = i
        n[1] = n[1]+1
    end
    
    sep_fit = Akima(angles[1,:], sep_vec[1,:])
    
    file2 = "../data/Hansen_Fig3_SepPoint.csv"
    hansen_sep_point = readdlm(file2, ',')
    
    error = zeros(1, 24)

    for i in 1:24
        error[1,i] = abs(hansen_sep_point[i+1, 2] - sep_fit(hansen_sep_point[i+1, 1].*pi/180))
    end

    average_error = sum(error)/length(error)

    @test average_error <= 0.01
end
