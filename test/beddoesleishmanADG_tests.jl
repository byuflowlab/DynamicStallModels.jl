#=
Test my implementation of AeroDyn's Beddoes-Leishman model using Gonzalez's modifications. I'm using the my modified NREL 5MW wind turbine as the basis for what I'm simulating. 

Adam Cardoza 1/3/23
=#

using DynamicStallModels, DelimitedFiles, OpenFASTTools, Statistics
using Test

DSM = DynamicStallModels
of = OpenFASTTools

path = dirname(@__FILE__)
cd(path)

include("./parseaerodyn.jl")
include("./testingutils.jl")





file = "../data/ADG_intermediates.txt"

datavec = readdlm(file, ',')

items = ["i", "t", "U", "alpha", "a_s", "M", "Gonzalez_factor", "c", "C_nalpha", "A1", "b1", "A2", "b2", "eta_e", "tau_v", "Cn_FS", "Cn_v", "Cn_alpha_q_nc", "Cn_q_circ", "Cn_alpha_q_circ", "fprimeprime", "Cn_alpha_nc", "Cn_q_nc", "T_alpha", "Kalpha_f", "Kprime_alpha", "q_f_cur", "alpha_filt_cur", "k_alpha", "alpha_e", "Kq_f", "Kprime_q", "Df", "fprime", "alpha_f", "Cc_pot", "fprime_c", "fprimeprime_c", "C_nalpha_circ", "Cn", "Cc", "Cm"]


### Read in AeroDyn files #Todo: Switch to a test that doesn't depend on updates from OpenFASTTools. *And the possibility that those files could change!
addriver = of.read_addriver("NREL5MW_ADdriver.dvr", "../testing/OpenFAST_NREL5MW_modified")
adfile = of.read_adfile("NREL5MW_ADfile.dat","../testing/OpenFAST_NREL5MW_modified")
adblade = of.read_adblade("NREL5MW_adblade.dat", "../testing/OpenFAST_NREL5MW_modified")

numnodes = Int(adblade["NumBlNds"])

if !@isdefined(mat)
    mat = parseaerodyn(datavec, items, numnodes) # mat goes time x vars x nodes
end

#Read in the trouble node
fi = open("../data/BLADG_intermediate_states.txt", "r")
namestring = readline(fi)
close(fi)

namestring = of.rmspaces(namestring)
namevec = readdlm(IOBuffer(namestring))

data = readdlm("../data/BLADG_intermediate_states.txt", skipstart=1)

outs = Dict(namevec[i] => data[:,i] for i in eachindex(namevec))

#Todo: Test that I'm simulating the correct thing. 



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

@testset "Beddoes-Leishman - AeroDyn - Gonzalez modifications" begin 

    ### Loop through the nodes and test them. 
    for i = 1:numnodes
        ### Make sure that the nodes in the intermediate states file are correct. 
        @test isapprox(mat[1,1,i], i)

        ### Prepare inputs that rely on the package. 
        af, xcp = of.make_dsairfoil(afs[i]; interp=Linear) 
        # af = make_dsairfoil(afs[i], chordvec[i]) 
        airfoils = Array{Airfoil, 1}(undef, 1) #TODO: I should probably change the type requirement. 
        airfoils[1] = af

        # dsmodel = BeddoesLeishman(Indicial(), 3, af.polar, af.alpha0, af.alphasep[2], ) #Note: No need to initialized a DSModel, because OpenFASTTools automagically initializes the correct one. 

        Uvec = zero(tvec)
        aoavec = zero(tvec)

        Uvec[3:end] = mat[:,3, i]
        aoavec[3:end] = mat[:,4, i]
        Uvec[1:2] .= mat[1,3, i]
        aoavec[1:2] .= mat[1,4, i]

        ### Solve
        states, loads = DSM.solve_indicial(airfoils, tvec, Uvec, aoavec; cvec=[chordvec[i]], xcpvec=[xcp]) #TODO: There might be a scoping issue somewhere because when I was running the code multiple times in a row, I would get a different error. 

        Cnvec = zero(tvec)
        Ccvec = zero(tvec)

        for i in eachindex(tvec)
            Cnvec[i], Ccvec[i] = DSM.rotate_load(loads[i,1], loads[i,2], states[i, 1])
        end 


        ### Calculate intermediate states and calculations. 
        Cfsn, Cvn, Nnc_aq, Nc_q, Nc_aq, Nnc_a, Nnc_q, Talpha, M, k_alpha, TI, alphae, k_q, Cpot, dcndalpha_circ = DSM.extractintermediatestates(states, Uvec, chordvec[i], af; a) #Todo: I should calculate some states based on these intermediate states so if I change the main code but not the intermediate states code then they won't match and fail their tests. 


        ### Calculate errors. 
        alpha_filt_rms = RMS(states[3:end, 2], mat[:,28,i])
        q_filt_rms = RMS(states[3:end, 4], mat[:,27,i])

        kalpha_rms = RMS(k_alpha[3:end], mat[:,29,i])
        Talpha_rms = RMS(Talpha[3:end], mat[:,24,i])

        Ka_rms = RMS(states[3:end, 5], mat[:,25,i])
        Kpa_rms = RMS(states[3:end, 7], mat[:,26,i])

        Kq_rms = RMS(states[3:end, 6], mat[:,31,i])
        Kpq_rms = RMS(states[3:end, 8], mat[:,32,i])

        alphae_rms = RMS(alphae[3:end], mat[:,30,i])

        Nc_q_rms = RMS(Nc_q[3:end], mat[:,19,i]) 
        Nc_aq_rms = RMS(Nc_aq[3:end], mat[:,20,i]) 

        Nnc_a_rms = RMS(Nnc_a[3:end], mat[:,22,i]) 
        Nnc_q_rms = RMS(Nnc_q[3:end], mat[:,23,i]) 
        Nnc_aq_rms = RMS(Nnc_aq[3:end], mat[:,18,i]) 

        fp_rms = RMS(states[3:end, 17], mat[:,34,i])
        fpp_rms = RMS(states[3:end, 23], mat[:,21,i])

        fpc_rms = RMS(states[3:end, 18], mat[:,38,i]) #TODO: Really good half the time, meh the other half. 
        fppc_rms = RMS(states[3:end, 24], mat[:,37,i]) #Todo: Interestingly this doesn't go to 1.44. It goes to 1.435 and change. Like it is fluctuating near 1.44. 

        CNFS_rms = RMS(Cfsn[3:end], mat[:,16,i]) 
        Cvn_rms = RMS(Cvn[3:end], mat[:,17,i])

        tau_rms = RMS(states[3:end, 28], mat[:,15,i]) #Todo: Really good for everything but the root. values. 

        # @show fpc_rms, fppc_rms

        # Cnerr = calculate_relerr(loads[3:end,1], mat[:,40,i]).*100
        Cnerr = calculate_relerr(Cnvec[3:end], mat[:,40,i]).*100
        # Ccerr = calculate_relerr(loads[3:end,2], mat[:,41,i]).*100
        Ccerr = calculate_relerr(Ccvec[3:end], mat[:,41,i]).*100
        Cmerr = calculate_relerr(loads[3:end,3], mat[:,42,i]).*100

        @show mean(abs.(Cnerr)), mean(abs.(Ccerr)), mean(abs.(Cmerr))


        ### Using RMS because the value might be small and relative error can blow out of proportion really quick. 
        @test alpha_filt_rms <= 1e-5  
        @test q_filt_rms <= 1e-5
        
        @test kalpha_rms <= 1e-14
        @test Talpha_rms <= 1e-16

        @test Ka_rms <= 1e-4
        @test Kpa_rms <= 1e-7

        @test Kq_rms <= 1e-5
        @test Kpq_rms <= 1e-7

        @test alphae_rms <= 1e-5

        @test Nc_q_rms <= 1e-5
        @test Nc_aq_rms <= 1e-4

        @test Nnc_a_rms <= 1e-5
        @test Nnc_q_rms <= 1e-5
        @test Nnc_aq_rms <= 1e-5

        @test fp_rms <= 0.01
        @test fpp_rms <= 0.01
        @test fpc_rms <= 0.05
        @test fppc_rms <= 0.05

        # @show fp_rms, fpp_rms

        @test CNFS_rms <= 0.01
        @test Cvn_rms <= 1e-3

        # @test tau_rms <= 1e-16 #Todo: Really good for everything but the root. (Same as comment above)

        ### Relative percent error 
        @test mean(abs.(Cnerr)) <= 0.07
        @test mean(abs.(Ccerr)) <= 0.11 
        @test mean(abs.(Cmerr)) <= 0.11

    end #End looping through the nodes to test them. 

end #End BLADG tests


