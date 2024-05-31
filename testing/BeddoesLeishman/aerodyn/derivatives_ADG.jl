using DynamicStallModels, OpenFASTTools, FLOWMath, Plots, Plots.PlotMeasures, DelimitedFiles, LaTeXStrings, Revise, Statistics

#=
Testing to see if I can't pass derivativess through the dynamic stall models. 

Adam Cardoza 7/17/23
=#



of = OpenFASTTools
DSM = DynamicStallModels

path = dirname(@__FILE__)
cd(path)

include("parseaerodyn.jl")

# items = ["t", "alpha", "Cc_pot", "fprime_c", "fprimeprime_c", "Df", "Cn", "Cc", "PI"] 
items = ["i", "t", "U", "alpha", "a_s", "M", "Gonzalez_factor", "c", "C_nalpha", "A1", "b1", "A2", "b2", "eta_e", "tau_v", "Cn_FS", "Cn_v", "Cn_alpha_q_nc", "Cn_q_circ", "Cn_alpha_q_circ", "fprimeprime", "Cn_alpha_nc", "Cn_q_nc", "T_alpha", "Kalpha_f", "Kprime_alpha", "q_f_cur", "alpha_filt_cur", "k_alpha", "alpha_e", "Kq_f", "Kprime_q", "Df", "fprime", "alpha_f", "Cc_pot", "fprime_c", "fprimeprime_c", "C_nalpha_circ", "Cn", "Cc", "Cm"] #, "Cn_visc", "Cc_visc"]
file = "/Users/adamcardoza/.julia/dev/DynamicStallModels/data/ADG_intermediates.txt"

entries = readdlm(file, ',')

mat = parseaerodyn(entries, items, 19)




fullout = readdlm("../../../data/gonzalez2.csv", ',')

names = fullout[1,:]
names[1] = "Time"
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in 1:length(names)) #11 or 12th node is what we're looking for. 

### Create airfoil
# c = 3.542
idx = 1



du21_a17 = of.read_airfoilinput("../../../data/airfoils/DU40_A17.dat")

dcndalpha = du21_a17.c_nalpha
alpha0 = du21_a17.alpha0*(pi/180)
polar = hcat(du21_a17.aoa.*(pi/180), du21_a17.cl, du21_a17.cd, du21_a17.cm)
clfit = Akima(polar[:,1], polar[:,2])
cdfit = Akima(polar[:,1], polar[:,3])
cmfit = Akima(polar[:,1], polar[:,4])

alphasep = [du21_a17.alpha2, du21_a17.alpha1].*(pi/180)

A = [du21_a17.a1, du21_a17.a2, du21_a17.a5]
b = [du21_a17.b1, du21_a17.b2]
T = [du21_a17.t_p, du21_a17.t_f0, du21_a17.t_v0, du21_a17.t_vl]

S1 = du21_a17.s1
S2 = du21_a17.s2
S3 = du21_a17.s3
S4 = du21_a17.s4
S = [S1, S2, S3, S4]
xcp = 0.2



airfoils = Array{Airfoil, 1}(undef, 1)





### Create model
# dsmodel = BeddoesLeishman(Indicial(), 1, airfoils, 3)


### Create simulation data
tvec = outs["Time"]
nt = length(tvec)

if idx <10
    num = "00$idx"
elseif idx<100
    num ="0$idx"
else
    num = "$idx"
end


Uvec = outs["AB1N"*num*"VRel"]
alphavec = (outs["AB1N"*num*"Alpha"].+0.01334).*(pi/180)

Uvec[3:end] = mat[:,3, idx]
alphavec[3:end] = mat[:,4, idx]

Udotvec = zero(Uvec)
alphadotvec = zero(alphavec)


### Solve

function objective(x)

    c , A1, A2, A5 = x

    af = of.make_dsairfoil(du21_a17, c; A=[A1, A2, A5])

    airfoils[1] = af

    states, loads = solve_indicial(airfoils, tvec[1:5], Uvec[1:5], alphavec[1:5])

    Cl = loads[:,1]

    return mean(Cl)
end

c0 = 3.542


x0 = [c0, A[1], A[2], A[3]]

Clavg = objective(x0)

using FiniteDiff, ForwardDiff

dCldx_fd = FiniteDiff.finite_difference_jacobian(objective, x0)
dCldx_fad = ForwardDiff.gradient(objective, x0)

@show dCldx_fd
@show dCldx_fad

nothing