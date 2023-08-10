using DynamicStallModels, DelimitedFiles, FLOWMath, DifferentialEquations, Test, Plots


dsm = DynamicStallModels


path = dirname(@__FILE__)
cd(path)

file = "../polars/Hansen_Figure3.csv"

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


parameters = [Uvector, 0.0, alphavec, alphavecdot]
x_initial = [0.0, 0.0, 0.0, 0.0]

prob = ODEProblem(airfoils, x_initial, tspan, parameters)
sol = DifferentialEquations.solve(prob, reltol=1e-8)

answer = parsesolution(dsmodel, airfoils, sol, parameters)


function hansen_fully_sep(airfoil , alpha)
    CL_st = airfoil.cl(alpha) #Finds the static lift coefficient values at the desired angle of attack
    dcldalpha = airfoil.dcldalpha
    alpha0 = airfoil.alpha0
    delta_alpha = (alpha-alpha0) #Used to multiply with the linear slope for inviscid lift

    fst = (2*sqrt(CL_st/(dcldalpha*(delta_alpha))) - 1)^2 #This equation is needed for the fully seperated lift equation to not give undefined results

    CL_fs = (CL_st - dcldalpha*delta_alpha*fst)/(1-fst)

    return CL_fs
end



file2 = "../polars/Hansen_Fig3_FullSep.csv"
hansen_separated = readdlm(file2, ',')
sep_fit = Akima(hansen_separated[:,1], hansen_separated[:,2])
error = zeros(1, 53)

@testset "Riso" begin
    @testset "Fully Separated Lift" begin
        for i in 0:52
            error[1, i+1] = abs(sep_fit((i/100)*180/pi) - hansen_fully_sep(af, i/100))
        end
        
        @test sum(error)/length(error) <= 0.02

        for i in 0:52
            @test hansen_fully_sep(af, i/100) <= af.cl(i/100)
        end
    end
end

