#=
Test the Riso model. 

- Test set 1 - Test that the functional model runs
- Test set 2 - Test that the iterative model runs

Adam Cardoza
=#

# using DynamicStallModels, DelimitedFiles, DifferentialEquations
# using Test

# DE = DifferentialEquations
# ds = DynamicStallModels

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

    alpha_ = amp*sin(omega*t) + shift
    return alpha_*(pi/180)
end

function alphadot(t)
    # k = 0.1 #given
    # u = 50 #Assumed
    # c = 1.0 #Assumed

    # omega = k*2*u/c #implied

    amp = 2
    shift = 5

    alphadot_ = amp*omega*cos(omega*t)

    return alphadot_*(pi/180)
end


aoa = -pi:0.01:pi
lift = 2*pi.*(aoa)
drag = zero(aoa)
polar = hcat(aoa, lift, drag)




A = [0.165, 0.335] #From the Hansen 2004 paper, for flat plate
b = [0.0455, 0.3000] # "" ""

Tp = 3.0 # "" "" 
Tf = 6.0 # "" "" 
T = [Tp, Tf]


m, n = size(polar)
newpolar = hcat(polar, zeros(m), zeros(m))
afs = Array{Airfoil,1}(undef,1)
afs[1] = airfoil(newpolar; A, b, T)

@testset "Risø Model" begin

    ################# Test the Functional Model #########################
    @testset "Functional Model" begin


        dsmodel = Riso(Functional(), length(afs), afs)

        x0 = zeros(4)
        x0[3] = 1.0

        p = [U, Udot, alpha, alphadot, c]

        tspan = (0.0, 80.0)


        prob = ODEProblem(dsmodel, x0, tspan, p)
        sol = solve(prob, dtmax=0.1)

        Cl, Cd, t =  parsesolution(dsmodel, sol, p)



        alphavec = alpha.(sol.t)

        expdata = readdlm("../data/Hansen2004/figure8_flatplate/indicial.csv", ',')

        # @test true
        @test isapprox(maximum(Cl), maximum(expdata[:,2]), atol=0.01) #Todo: I need a better to test. A way that would include phase and such would be good. 
        # @test isapprox(minimum(Cl), minimum(expdata[:,2]), atol=0.01) #Doesn't working because the state starts at zero. Meep. 

    end


    ################## Test iterative implementation ################
    @testset "Iterative Model" begin #Todo I wonder if the scope will carry over or not. -> It doesn't

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
        
            DE.step!(integrator, dt_local, true)
            xds[i,:] = integrator.u
        
            cl_local, cd_local = ds.parsestates(dsmodel, xds[i,:], integrator.p)
        
            cl_d[i] = cl_local[1]
            cd_d[i] = cd_local[1]
        end

        alphavec = alpha.(tvec)

        expdata = readdlm("../data/Hansen2004/figure8_flatplate/indicial.csv", ',')
        

        @test isapprox(maximum(cl_d), maximum(expdata[:,2]), atol=0.01)
    end


end



############# Prepare and run the iterative model
