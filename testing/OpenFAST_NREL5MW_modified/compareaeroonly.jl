import OpenFASTsr
using Plots, FLOWMath, Statistics, DelimitedFiles, CCBlade

of = OpenFASTsr
fm = FLOWMath

path = dirname(@__FILE__)
cd(path)

filename = "NREL5MW_aeroonly.1.out"



fullout = readdlm(filename; skipstart=6)

names = fullout[1,:]
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in 1:length(names))

adblade = of.read_adblade("NREL5MW_adblade.dat", path)
adfile = of.read_adfile("NREL5MW_ADfile.dat", path)
addriver = of.read_addriver("NREL5MW_ADdriver.inp", path)


rvec = adblade.span .+ addriver.hubrad
chordvec = adblade.chord
twistvec = adblade.twist .* (pi/180)


nt = length(outs["Time"])
na = length(rvec)


fxmat = zeros(nt, na)
fymat = zeros(nt, na)

for i = 1:na
    if i<10
        number = "00$i"
    elseif i<100
        number = "0$i"
    else
        number = "$i"
    end
    namex = "AB1N"*number*"Fx"
    namey = "AB1N"*number*"Fy"
    fxmat[:,i] = outs[namex]
    fymat[:,i] = outs[namey]
end


tvec = outs["Time"]

rhub = addriver.hubrad
rtip = adblade.span[end] + rhub
rotor = Rotor(rhub, rtip, 3, precone=0.0, turbine=true)

aftypes = Array{AlphaAF}(undef, 8)
aftypes[1] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder1.dat", radians=false)
aftypes[2] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder2.dat", radians=false)
aftypes[3] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU40_A17.dat", radians=false)
aftypes[4] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU35_A17.dat", radians=false)
aftypes[5] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU30_A17.dat", radians=false)
aftypes[6] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU25_A17.dat", radians=false)
aftypes[7] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU21_A17.dat", radians=false)
aftypes[8] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/NACA64_A17.dat", radians=false)

airfoils = aftypes[adblade.afid]

sections = Section.(rvec, chordvec, twistvec, airfoils)

yaw = addriver.yaw[1]
tilt = addriver.shfttilt[1]
hubHt = addriver.hubht
shearExp = addriver.shearexp[1]

Vinf = addriver.windspeed[1]
rotorR = rtip*cos(precone)

Omega = addriver.rpm[1]*2*pi/60
pitch = addriver.pitch[1]
azimuth = 0.0*pi/180
rho = 1.225

op = windturbine_op.(Vinf, Omega, pitch, rvec, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)

out = solve.(Ref(rotor), sections, op)



forceplt = plot(xaxis="Radius (m)", yaxis="Distributed Load (N/m)", legend=:topleft)
plot!(rvec, fxmat[end,:], lab="Fx")
plot!(rvec, fymat[end,:], lab="Fy")
plot!(rvec, out.Np, lab="CCBlade N" )
plot!(rvec, out.Tp, lab="CCBlade T")
display(forceplt)

# writedlm("fxmat.csv", fxmat, ',' )
# writedlm("fymat.csv", fymat, ',' )

#=
Adam Cardoza 7/27/22
Well they match decently well. I think the main difference is that AeroDyn fluctuates with time and they likely use slightly different polars... which would explain the differences. I suppose the assumed density could be different, but I don't think that woould make a huge difference. 

Interesting that the loading is like 1000 Newtons less at the peak than the static solution that includes fewer points... Like... it's the exact same airfoils... 

Oh... the blade is different. 
=#

