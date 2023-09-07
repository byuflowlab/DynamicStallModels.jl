using DynamicStallModels, DelimitedFiles, DifferentialEquations, OpenFASTTools, FLOWMath, Statistics
using Test

dsm = DynamicStallModels
of = OpenFASTTools

path = dirname(@__FILE__)
cd(path)

include("./testingutils.jl")


@testset "Airfoils" begin

    @testset "Airfoil Constructors" begin
        @test true
    end #End Airfoil Constructor tests

    @testset "Airfoil Methods" begin
        let  ### Test Airfoil ability to update states. 
            c = 0.1

            M = 0.379
            a = 343.0
            Vrel = M*a 

            dsmodel = Oye(Indicial(), 1, 1, 4.0)

            du21_a17 = of.read_airfoilinput("../data/airfoils/DU40_A17.dat") 
            af = of.make_dsairfoil(du21_a17, c) 
            af = update_airfoil(af; dsmodel, sfun=dsm.OSP(), alphasep=af.alphasep.*4.5)

            airfoils = Array{Airfoil, 1}(undef, 1)
            airfoils[1] = af


            tvec = range(0, 0.05, 1000) #0:0.001:0.05
            Uvec = Vrel.*ones(length(tvec))

            function alpha(t)
                c = 0.1
                M = 0.379
                a = 343.0
                shift = 10.3
                amp = 8.1
                k = 0.075
            
                v = M*a
                omega = k*2*v/c
            
                alf = shift + amp*sin(omega*t)
                return alf*(pi/180)
            end

            alphavec = alpha.(tvec)

            states, loads = solve_indicial(airfoils, tvec, Uvec, alphavec) #By running solve_indicial, I test   the airfoil method to call the update_states function of the aforementioned dynamic stall model. 

            @test mean(states)!=0.0
            @test !any(i -> isnan(i), states) #Test for NaN
            @test !any(i -> isinf(i), states) #Test for Inf
        end

        let ### Test vectors of airfoils functionalities
            ## Read in AeroDyn files #Todo: Switch to a test that doesn't depend on updates from OpenFASTTools.
            addriver = of.read_addriver("NREL5MW_ADdriver.dvr", "../testing/OpenFAST_NREL5MW_modified")
            adfile = of.read_adfile("NREL5MW_ADfile.dat","../testing/OpenFAST_NREL5MW_modified")
            adblade = of.read_adblade("NREL5MW_adblade.dat", "../testing/OpenFAST_NREL5MW_modified")

            numnodes = Int(adblade["NumBlNds"])

            #### Define variables. 
            indices = 1:numnodes
            chordvec = adblade["BlChord"][indices] 
            twistvec = (pi/180) .* adblade["BlTwist"][indices] 
            rho = addriver["FldDens"] 
            mu = addriver["KinVisc"] 
            a = addriver["SpdSound"]

            U = 100.0 #m/s
            Udot = 0.0 
            alpha_1 = 2*pi/180
            alpha_2 = 2*pi/180
            alphadot = 0.0

            ### Prep the airfoils 
            aftypes = Array{of.AirfoilInput}(undef, 8)
            aftypes[1] = of.read_airfoilinput("../testing/OpenFAST_NREL5MW_modified/Airfoils/Cylinder1.dat") 
            aftypes[2] = of.read_airfoilinput("../testing/OpenFAST_NREL5MW_modified/Airfoils/Cylinder2.dat") 
            aftypes[3] = of.read_airfoilinput("../testing/OpenFAST_NREL5MW_modified/Airfoils/DU40_A17.dat") 
            aftypes[4] = of.read_airfoilinput("../testing/OpenFAST_NREL5MW_modified/Airfoils/DU35_A17.dat") 
            aftypes[5] = of.read_airfoilinput("../testing/OpenFAST_NREL5MW_modified/Airfoils/DU30_A17.dat") 
            aftypes[6] = of.read_airfoilinput("../testing/OpenFAST_NREL5MW_modified/Airfoils/DU25_A17.dat") 
            aftypes[7] = of.read_airfoilinput("../testing/OpenFAST_NREL5MW_modified/Airfoils/DU21_A17.dat") 
            aftypes[8] = of.read_airfoilinput("../testing/OpenFAST_NREL5MW_modified/Airfoils/NACA64_A17.dat")


            # indices correspond to which airfoil is used at which station
            af_idx = Int.(adblade["BlAFID"][indices]) 


            # create airfoil array
            afs = aftypes[af_idx]

            dt = addriver["DT"][1] 
            tvec = [0.0, 0.1]

            airfoils = [of.make_dsairfoil(afs[i], chordvec[i]; interp=Linear) for i in eachindex(afs)]

            n = length(airfoils)

            ### Initialize the data storage
            ns = dsm.numberofstates_total(airfoils)

            states = Array{eltype(U), 2}(undef, 2, ns)
            y = Array{eltype(U), 1}(undef, 4n)
            for k = 1:n
                kdxs = 4*(k-1)+1:4k
                y[kdxs] = [U, Udot, alpha_1+twistvec[k], alphadot]
            end
            stateidx = Vector{Int}(undef, n)

            tempx = 1
            for i in eachindex(airfoils)
                stateidx[i] = tempx
                nsi1, nsi2 = dsm.state_indices(airfoils[i].model, stateidx[i])
                loadidx = 3*(i-1)+1:3i
                paramidx = 4*(i-1)+1:4*i 
                ys = view(y, paramidx)
                states[1,nsi1:nsi2], _ = dsm.initialize(airfoils[i], tvec, ys) 
                tempx += dsm.numberofstates(airfoils[i].model)
            end

            for k = 1:n
                kdxs = 4*(k-1)+1:4k
                y[kdxs] = [U, Udot, alpha_2+twistvec[k], alphadot]
            end

            # @show y

            xsi = view(states, 1, :) #Note: You can't pass in a slice, you must pass in a view to this function. 
            xs1 = view(states, 2, :)
            airfoils(xsi, xs1, stateidx, y, dt)


            @test !any(item -> isnan(item), states)
            @test mean(states[2,:])!=0.0
            @test !any(i -> isinf(i), states) #Test for Inf


        end # End testing vectors of airfoils


    end #End Airfoil Method tests

    @testset "Separation Point Functions" begin
        @test true
    end #End Airfoil separation point function tests

end