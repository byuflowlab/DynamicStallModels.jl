#=
Test my implementation of AeroDyn's Beddoes-Leishman model using Gonzalez's modifications. 

Adam Cardoza 1/3/23
=#

using DynamicStallModels, DelimitedFiles, OpenFASTsr, Statistics
using Test

DSM = DynamicStallModels
of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)

include("./parseaerodyn.jl")

err(x, xt) = x-xt
relerr(x, xt) = (x-xt)/xt
function RMS(x, xt)
    diff = @. (x -xt)^2
    return sqrt(sum(diff)/length(x))
end




file = "../data/ADG_intermediates.txt"

datavec = readdlm(file, ',')

items = ["i", "t", "U", "alpha", "a_s", "M", "Gonzalez_factor", "c", "C_nalpha", "A1", "b1", "A2", "b2", "eta_e", "tau_v", "Cn_FS", "Cn_v", "Cn_alpha_q_nc", "Cn_q_circ", "Cn_alpha_q_circ", "fprimeprime", "Cn_alpha_nc", "Cn_q_nc", "T_alpha", "Kalpha_f", "Kprime_alpha", "q_f_cur", "alpha_filt_cur", "k_alpha", "alpha_e", "Kq_f", "Kprime_q", "Df", "fprime", "alpha_f", "Cc_pot", "fprime_c", "fprimeprime_c", "C_nalpha_circ", "Cn", "Cc"]


### Read in AeroDyn files
addriver = of.read_addriver("NREL5MW_ADdriver.dvr", "../testing/OpenFAST_NREL5MW_modified")
adfile = of.read_adfile("NREL5MW_ADfile.dat","../testing/OpenFAST_NREL5MW_modified")
adblade = of.read_adblade("NREL5MW_adblade.dat", "../testing/OpenFAST_NREL5MW_modified")

numnodes = Int(adblade["NumBlNds"])

if !@isdefined(mat)
    mat = parseaerodyn(datavec, items, numnodes) # mat goes time x vars x nodes
end



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

tspan = (0.0, addriver["Tmax"][1]) 
dt = addriver["dT"][1] 
tvec = tspan[1]:dt:4.9

@testset "Beddoes-Leishman - AeroDyn - Gonzalez modifications" begin 

    ### Loop through the nodes and test them. 
    for i = 1:numnodes

        @test isapprox(mat[1,1,i], i)

        af, constants = make_dsairfoil(afs[i])
        airfoils = Array{Airfoil, 1}(undef, 1) #TODO: I should probably change the type requirement. 
        airfoils[1] = af

        dsmodel = BeddoesLeishman(Indicial(), 1, airfoils, 3, constants)

        Uvec = zero(tvec)
        aoavec = zero(tvec)

        Uvec[3:end] = mat[:,3, i]
        aoavec[3:end] = mat[:,4, i]
        Uvec[1:2] .= mat[1,3, i]
        aoavec[1:2] .= mat[1,4, i]

        ### Solve
        states, loads = solve_indicial(dsmodel, [chordvec[i]], tvec, Uvec, aoavec; a)

        alpha_filt_rms = RMS(states[3:end,1], mat[:,28,i])
        q_filt_rms = RMS(states[3:end,29], mat[:,27,i])

        Cnerr = relerr(loads[3:end,1], mat[:,40,i]).*100
        Ccerr = relerr(loads[3:end,2], mat[:,41,i]).*100

        # @show mean(abs.(Cnerr)), mean(abs.(Ccerr))

        ### Using RMS because the value might be small and relative error can blow out of proportion really quick. 
        @test alpha_filt_rms <= 1e-5  
        @test q_filt_rms <= 1e-5  

        ### Going back to relative percent error
        @test mean(abs.(Cnerr)) <= 0.1
        @test mean(abs.(Ccerr)) <= 0.1

    end #End looping through the nodes to test them. 

end #End BLADG tests


