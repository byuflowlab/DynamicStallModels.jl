using DynamicStallModels, DelimitedFiles, DifferentialEquations, OpenFASTsr, FLOWMath, Statistics
using Test

#Todo: This appears to be broken. 

DE = DifferentialEquations
ds = DynamicStallModels
of = OpenFASTsr

fullout = readdlm("../data/aerodynout_fordynamicstall.csv", ',')
names = fullout[1,:]
names[1] = "Time"
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in 1:length(names)) #11 or 12th node is what we're looking for. 

### Create airfoil
c = 1.0

du21_a17 = of.read_airfoilinput("../data/DU21_A17.dat")

dcndalpha = du21_a17.c_nalpha
alpha0 = du21_a17.alpha0*(pi/180)
polar = hcat(du21_a17.aoa.*(pi/180), du21_a17.cl, du21_a17.cd, du21_a17.cm)
clfit = Akima(polar[:,1], polar[:,2])
cdfit = Akima(polar[:,1], polar[:,3])
cmfit = Akima(polar[:,1], polar[:,4])

alphasep = [du21_a17.alpha2, du21_a17.alpha1].*(pi/180)
# alphasep = [-8, 14].*(pi/180) #Flow was seperating when I didn't expect it to. 

A = [du21_a17.a1, du21_a17.a2]
b = [du21_a17.b1, du21_a17.b2]
T = [du21_a17.t_p, du21_a17.t_f0, du21_a17.t_v0, du21_a17.t_vl]

S1 = du21_a17.s1
S2 = du21_a17.s2
S3 = du21_a17.s3
S4 = du21_a17.s4
S = [S1, S2, S3, S4]
xcp = 0.2

af = Airfoil(polar, clfit, cdfit, cmfit, dcndalpha, alpha0, alphasep, A, b, T, S, xcp)
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af

### Create simulation data
tvec = outs["Time"]
nt = length(tvec)

Uvec = [sqrt(outs["AB1N011Vx"][i]^2 + outs["AB1N011Vy"][i]^2) for i in 1:nt] #m/s
aoavec = outs["AB1N011Alpha"].*(pi/180)

@testset "Beddoes-Leishman Model" begin
    @testset "AeroDyn Original Implementation" begin


        ### Create model
        zeta = -3 #Low-pass-filter frequency cutoff
        A5 = du21_a17.a5
        b5 = du21_a17.b5
        Tsh = du21_a17.st_sh #Strouhal Frequency
        eta = du21_a17.eta_e #Recovery factor
        constants = [zeta, A5, b5, Tsh, eta]
        dsmodel = BeddoesLeishman(Indicial(), 1, airfoils, 2, constants)


        ### Solve
        states, loads = solve_indicial(dsmodel, [c], tvec, Uvec, aoavec)

        Cn = loads[:,1]

        @test isapprox(maximum(outs["AB1N011Cn"]), maximum(Cn), atol=0.01)
        @test isapprox(minimum(outs["AB1N011Cn"]), minimum(Cn), atol=0.01)
        @test isapprox(mean(outs["AB1N011Cn"]), mean(Cn), atol=0.01)
    end
end
