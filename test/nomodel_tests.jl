#=
Testing NoModel. 
=#

using Test
using DynamicStallModels, DelimitedFiles, OpenFASTsr, Statistics

DSM = DynamicStallModels
of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)

include("./testingutils.jl")


@testset "NoModel" begin
    ### Read in AeroDyn files #Todo: Switch to a test that doesn't depend on updates from OpenFASTsr.
    #Todo: I might have changed the input files... which would be a bummer. 
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


    ### Prep the ASD rotor and operating conditions 
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

    tspan = (0.0, addriver["TMax"][1]) 
    dt = addriver["DT"][1] 
    tvec = tspan[1]:dt:4.9

    ### Prepare inputs that rely on the package. 
    
    airfoils = Array{Airfoil, 1}(undef, numnodes) 
    
    for i in eachindex(chordvec)
        airfoils[i] = make_dsairfoil(afs[i], chordvec[i]; interp=Linear) 
        if i<5
            airfoils[i] = update_airfoil(airfoils[i]; dsmodel=DSM.NoModel())
        end
    end

    #Check that we got mixed models in the airfoils. 
    for ti in eachindex(airfoils)
        if ti<5
            @test isa(airfoils[ti].model, DSM.NoModel)
        else
            @test isa(airfoils[ti].model, BeddoesLeishman)
        end
    end

    Vrel = 130.0

    nt = 100
    tvec = collect(range(0, 1, nt))
    Uvec = ones(nt).*Vrel

    function alpha(t)
        c = 0.1
        shift = 10.3
        amp = 8.1
        k = 0.075
    
        v = Vrel
        omega = k*2*v/c
    
        alf = shift + amp*sin(omega*t)
        return alf*(pi/180)
    end

    alphavec = alpha.(tvec)

    states, loads = solve_indicial(airfoils, tvec, Uvec, alphavec)


    for j in eachindex(airfoils)
        loadidx = 3*(j-1)+1:3j
        if j<5
            clvec = loads[:, loadidx[1]]
            cdvec = loads[:, loadidx[2]]
            cmvec = loads[:, loadidx[3]]

            clvec_gold = airfoils[j].cl.(alphavec)
            cdvec_gold = airfoils[j].cd.(alphavec)
            cmvec_gold = airfoils[j].cm.(alphavec)

            @test isapprox(clvec, clvec_gold)
            @test isapprox(cdvec, cdvec_gold)
            @test isapprox(cmvec, cmvec_gold)
        else
            @test mean(loads[:,loadidx])!=0.0
            #TODO: Should probably check that the no model doesn't affect the other states... or loads. 
        end
    end

end #End no model tests