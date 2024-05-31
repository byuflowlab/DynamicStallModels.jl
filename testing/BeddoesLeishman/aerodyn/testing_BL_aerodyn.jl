using DynamicStallModels, OpenFASTTools, FLOWMath, Plots, Plots.PlotMeasures

#=
Test the Beddoes-Leishman model as given in the documentation. I had to make some minor modifications as not all of the equations were given. 

Adam Cardoza 8/22/22
=#

of = OpenFASTTools

path = dirname(@__FILE__)
cd(path)


include("../../src/BeddoesLeishman/BeddoesLeishmanAeroDyn.jl")



### Create airfoil
c = 1.0

du21_a17 = of.read_airfoilinput("../../data/DU21_A17.dat")

dcndalpha = du21_a17.c_nalpha
alpha0 = du21_a17.alpha0*(pi/180)
polar = hcat(du21_a17.aoa.*(pi/180), du21_a17.cl, du21_a17.cd, du21_a17.cm)
clfit = Akima(polar[:,1], polar[:,2])
cdfit = Akima(polar[:,1], polar[:,3])
cmfit = Akima(polar[:,1], polar[:,4])

# alphasep = [du21_a17.alpha2, du21_a17.alpha1].*(pi/180)
alphasep = [-8, 14].*(pi/180) #Flow was seperating when I didn't expect it to. 

A = [du21_a17.a1, du21_a17.a2]
b = [du21_a17.b1, du21_a17.b2]
T = [du21_a17.t_p, du21_a17.t_f0, du21_a17.t_v0, du21_a17.t_vl]

af = Airfoil(polar, clfit, cdfit, cmfit, dcndalpha, alpha0, alphasep, A, b, T)
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af

S1 = du21_a17.s1
S2 = du21_a17.s2
S3 = du21_a17.s3
S4 = du21_a17.s4


### Create model
zeta = -3 #Low-pass-filter frequency cutoff
A5 = du21_a17.a5
b5 = du21_a17.b5
Tsh = du21_a17.st_sh #Strouhal Frequency
eta = du21_a17.eta_e #Recovery factor
constants = [zeta, A5, b5, Tsh, eta]
dsmodel = BeddoesLeishman(Indicial(), 1, airfoils, 1, constants)


### Create simulation data
U = 10 #m/s
tspan = (0,3)
dt = 0.1 #The model get's closer when increasing time step... that can't be right. 
aoa = 5*(pi/180) #degrees
amp = 3*(pi/180)
frequency = 1 #pi #2*pi #Decreasing the frequency decreases the maximum loading.... which is very odd. 

tvec = tspan[1]:dt:tspan[2]
nt = length(tvec)
Uvec = ones(nt).*U
aoavec = zero(tvec)

for i = 1:nt
    aoavec[i] = amp*sin(frequency*tvec[i]) + aoa
end





### Solve
states, Cn, Cc, Cl, Cd, Cm = solve_indicial(dsmodel, c, tvec, Uvec, aoavec, S1, S2, S3, S4; a = 343)


#=
If the states are changing, they don't appear to be changing enough to make a difference. But at least it runs! -> The final value is pretty close to what the static value is... so that's good I guess. It is like 8.9% off... which isn't great... but it is pretty close.

With an oscillating aoa, the bad boy explodes. -> I was initializing at a bad aoadot. 

It is interesting that the model's Cn leads the angle of attack. Interesting as in wrong. 

Well... In the very least, I could take the angle of attack output and Cn values from an Aerodyn run and input them into the function to see if the model is anywhere close to AeroDyn's dynamic stall model. 

=#

staticCn = clfit.(aoavec)

cnplt = plot(xaxis="time (s)", yaxis="Cn", right_margin=20mm)
plot!(tvec, Cn, lab="BL_AD", markershape=:x)
plot!(tvec, staticCn, lab="static")
plot!(twinx(), tvec, aoavec.*(180/pi), lab="AOA", leg=:bottomleft, linecolor=:green, linestyle=:dash, yaxis="AOA")
display(cnplt)

nothing