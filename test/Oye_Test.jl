using DynamicStallModels, DelimitedFiles, FLOWMath, DifferentialEquations, Test

dsm = DynamicStallModels

path = dirname(@__FILE__)
cd(path)

c = 0.55

M = 0.11
a = 343.0
Vrel = M*a 

file = "../polars/NACA_0015_Faber.csv"

file_2 = "../polars/Faber_0015_Oye_Top_Half.csv"

file_3 = "../polars/Faber_0015_Oye_Bottom_Half.csv"

polar = readdlm(file , ',') #reads in Faber's static lift polar for the NACA_0015 airfoil

Polar_Top_Half = readdlm(file_2, ',') #reads in the dynamic lift values from Faber's NACA 0015 (top half)
Lift_Top_Half = Akima(Polar_Top_Half[:,1].*pi/180, Polar_Top_Half[:,2])

Polar_Bottom_Half = readdlm(file_3, ',') #reads in the dynamic lift values from Faber's NACA 0015 (bottom half)
Lift_Bottom_Half = Akima(Polar_Bottom_Half[:,1].*pi/180, Polar_Bottom_Half[:,2])

dsmodel = Oye(Functional(), 1, 2, 4.0)

af = dsm.make_airfoil(polar, dsmodel, c; sfun=dsm.OSP())

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af

tspan = (0, 2.0) #0:0.001:0.05

function Uvector(t)
    return 0.11*343.0
end

function alpha(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.0
    k = 0.051

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

parameters = [Uvector, 0.0, alpha, 0.0]

x_initial = [0.8]

prob = ODEProblem(airfoils, x_initial, tspan, parameters)

sol = DifferentialEquations.solve(prob, reltol=1e-8)

answer = parsesolution(dsmodel, airfoils, sol, parameters)

dsmodel_2 = Oye(Indicial(), 1, 2, 4.0)

af_2 = dsm.make_airfoil(polar, dsmodel_2, c; sfun=dsm.OSP())

airfoils_2 = Array{Airfoil, 1}(undef, 1)
airfoils_2[1] = af_2

tvec = range(0, 1.0, 1000)
Uvec = Vrel.*ones(length(tvec))
alphavec = alpha.(tvec)


states, loads = solve_indicial(airfoils_2, tvec, Uvec, alphavec)

@testset "Oye" begin
    @testset "Oye Functional Test" begin
        #tests to see if the states are always between 0 and 1
        @testset "State Value" begin
            f = Array(sol)
            for i in 1:length(f)
                @test 0 <= f[i] <= 1
            end
        end
        #test to see if the fully separated lift is always positive and below the static lift; test to see if the derivative near 32 degrees is pi/6
        #test to see if the derivative near alpha0 is pi
        @testset "Fully Separated Lift" begin
            for i in 1:length(sol.t)
                C_fs = dsm.cl_fullysep_faber(af, alpha(sol.t[i]))
                Cl = af.cl(alpha(sol.t[i]))
                @test C_fs <= Cl
                @test C_fs >= -0.05 #since they might get slightly negative close to an aoa of zero
            end


            C_fs_end = dsm.cl_fullysep_faber(af, 32*pi/180)
            C_fs_end_minus_step = dsm.cl_fullysep_faber(af, 31.99*pi/180)


            derivative_1 = (C_fs_end-C_fs_end_minus_step)/(32*pi/180 - 31.99*pi/180) 


            @test isapprox(derivative_1, (pi)/6 , rtol = 0.05) #testing derivative at alphasep


            C_fs_beginning = dsm.cl_fullysep_faber(af, af.alpha0)
            C_fs_beginning_plus_step = dsm.cl_fullysep_faber(af, af.alpha0+0.01)


            derivative_2 = (C_fs_beginning_plus_step- C_fs_beginning)/(0.01)


            @test isapprox(derivative_2, pi, rtol=0.05) #testing derivative at alpha0
        end
        @testset "Dynamic Lift Compared to Faber (Functional)" begin
           for i in 1:length(answer[1,1:(end-1)]) #compares the lfit from Faber to our Lift using the Functional Method method
                if answer[1,i] < answer[1, i+1]
                    Faber_Lift = Lift_Top_Half(answer[1,i])
                else
                    Faber_Lift = Lift_Bottom_Half(answer[1,i])
                end
                @test isapprox(answer[2,i], Faber_Lift, atol = 0.05)
            end
        end
     end

   
    @testset "Oye Indicial Test" begin
        #test to see if the state value is between 0 and 1
        @testset "State Values" begin
            for i in 1:length(states)
                @test 0 <= states[i] <= 1
            end
        end

        @testset "Fully Separated Lift" begin
            for w in 1:length(airfoils_2)
                for i in 1:length(tvec)
                    Cl_fs = dsm.cl_fullysep_faber(airfoils_2[w], alphavec[i])
                    Cl = airfoils_2[w].cl(alphavec[i])
                    @test Cl_fs <= Cl #tests to see if the fully sep lift is always less than or equal to the static lift
                    @test Cl_fs >= -0.05 #tests to see if fully sep lift is greater than zero (near alpha0 it might be a little bit negative)
                end


                C_fs_end = dsm.cl_fullysep_faber(airfoils_2[w], 32*pi/180)
                C_fs_end_minus_step = dsm.cl_fullysep_faber(airfoils_2[w], 31.99*pi/180)


                derivative_1 = (C_fs_end-C_fs_end_minus_step)/(32*pi/180 - 31.99*pi/180) 


                @test isapprox(derivative_1, (pi)/6 , rtol = 0.05) #testing derivative at alphasep


                C_fs_beginning = dsm.cl_fullysep_faber(airfoils_2[w], af.alpha0)
                C_fs_beginning_plus_step = dsm.cl_fullysep_faber(airfoils_2[w], af.alpha0+0.01)


                derivative_2 = (C_fs_beginning_plus_step- C_fs_beginning)/(0.01)


                @test isapprox(derivative_2, pi, rtol=0.05) #testing derivative at alpha0
            end
        end

        @testset "Dynamic Lift Compared to Faber (Indicial)" begin
            for i in 1:(length(tvec)-1) #compares the lfit from Faber to our Lift using the Indicial method
                if alphavec[i] < alphavec[i+1]
                    Faber_Lift = Lift_Top_Half(alphavec[i])
                else
                    Faber_Lift = Lift_Bottom_Half(alphavec[i])
                end
                @test isapprox(loads[i,1], Faber_Lift, atol = 0.05)
            end
        end
    end
end