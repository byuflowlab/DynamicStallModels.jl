import OpenFASTsr
using Plots, FLOWMath, Statistics, DelimitedFiles, CCBlade

#=
Adam Cardoza 7/25/22

Prepare files for an OpenFAST simulation of a modified version of the NREL 5MW. No precone, no tilt, no yaw. Assuming a blade with a thickness 8% of the chord length. The beam will be modeled as a rectangular aluminum beam. 

=#

function getfieldnames(obj)
    return fieldnames(typeof(obj))
end

of = OpenFASTsr
fm = FLOWMath

path = dirname(@__FILE__)
cd(path)

#### Define variables. 
## Define simplified NREL 5MW Turbine constants and other info. 

density=3.0e2 #2.69e3
E= 7.0e8 #6.83e10
nu=0.34
G = E/(2*(1+nu))

thickness = 0.08
rhub = 1.5
rtip = 63.0
rvec = [11.7500, 15.8500, 19.9500, 24.0500, 28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333]
chordvec = [4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
twistvec = pi/180*[13.308, 11.480, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
thickvec = chordvec.*thickness
B = 1.0
hubht = 90.0


## Airfoil Constants
A = [0.3, 0.7] #Aerodyn defaults
b = [0.13, 0.53]
Tp = 1.7
Tf = 3.0


## Turbine Control variables
pitch = 0.0
precone = 0.0 #2.5*pi/180 #Todo: !!!! I need to work in a way to include precone
yaw = 0.0*pi/180
tilt = 0.0 #5.0*pi/180
azimuth = 0.0

## Environmental variables
vinf = 10.0
tsr = 7.55
rotorR = rtip*cos(precone)
omega = vinf*tsr/rotorR
rpm = omega*60/(2*pi)


frequency = 1.0
amplitude = 0.0
rho = 1.225
mu = 18.13e-6
a = 343.0
shearexp = 0.0

#### Define solution
tspan = (0.0, 0.5)
dt = 0.01

### Prep the ASD rotor and operating conditions 
aftypes = Array{AlphaAF}(undef, 8)
aftypes[1] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder1.dat", radians=false)
aftypes[2] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/Cylinder2.dat", radians=false)
aftypes[3] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU40_A17.dat", radians=false)
aftypes[4] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU35_A17.dat", radians=false)
aftypes[5] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU30_A17.dat", radians=false)
aftypes[6] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU25_A17.dat", radians=false)
aftypes[7] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/DU21_A17.dat", radians=false)
aftypes[8] = AlphaAF("/Users/adamcardoza/.julia/dev/CCBlade/data/NACA64_A17.dat", radians=false)

# indices correspond to which airfoil is used at which station
af_idx = [3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]

# create airfoil array
airfoils = aftypes[af_idx]

# p_s, xp, xe = create_simplebeam(rvec, chordvec, twistvec, rhub, rtip, thickvec)

# ## Create models
# gxmodel = gxbeam(xp, xe)







###### Read in files from disk
adfile = of.read_adfile("NREL5MW_ADfile.dat", path)
adblade = of.read_adblade("NREL5MW_adblade.dat", path)

edfile = of.read_edfile("NREL5MW_EDfile.dat", path) 
edblade = of.read_edblade("NREL5MW_EDblade.dat", path)

bdfile = of.read_bdfile("NREL5MW_BDfile.dat", path)
bdblade = of.read_bdblade("NREL5MW_BDblade.dat", path)



######## Update ADFiles
adfile.twraero = false
adfile.ai_drag = true
adfile.ti_drag = true
adfile.blades[:] .= "NREL5MW_adblade.dat"

adrads = adblade.span .+ edfile.hubrad


######### Update Elastodyn files
nedr = length(edblade.frac)

edfile.rotspeed = rpm
edfile.blpitchs[:] .= pitch
edfile.precones[:] .= precone
edfile.bldnodes = nedr
edfile.bldfiles[:] .= "NREL5MW_EDblade.dat"
edfile.twrfile = "NREL5MW_EDtower.dat"




edblade.adjblms = 1.057344 #AD15 factor

edradius = edblade.frac.*(edfile.tiprad - edfile.hubrad) .+ edfile.hubrad
chordfit = Akima(adrads, adblade.chord)
edchords = chordfit.(edradius)
edthicks = edchords.*thickness

Iz = zeros(nedr)
Iy = zeros(nedr)
rho = zeros(nedr) #Distributed weight. 
for i = 1:nedr
    w = edchords[i] 
    h = edthicks[i] 
    Area = h*w
    rho[i] = density*Area
    Iz[i] = E*w*(h^3)/12 #Second moment of area
    Iy[i] = E*h*(w^3)/12 #TODO: blade mode shapes will be off. 
end
edblade.flapstiff = Iz #TODO: I'm not 100% sure of this. I really neeed to look it up and figure it out. 
edblade.edgestiff = Iy


######### Update BDblade
nb = length(bdblade.nodes)

bdfracs = [bdblade.nodes[i].frac for i=1:nb]
bdradius = bdfracs.*(edfile.tiprad - edfile.hubrad) .+ edfile.hubrad
bdchords = chordfit.(edradius)
bdthicks = edchords.*thickness



for i = 1:nb
    massmat = zeros(6,6)
    stiffmat = zeros(6,6)

    w = bdchords[i] 
    h = bdthicks[i] 
    Area = h*w
    rho_local = density*Area
    massmat[1,1] = massmat[2,2] = massmat[3,3] = rho_local #TODO: will all of massmat change w/ rho_local? Like.... if I change rho_local on iteration 2, will what I store massmat on iteration 1 also change? 
    Iz_local = w*(h^3)/12 #Second moment of area - Flatwise I believe
    Iy_local = h*(w^3)/12 # SMOA - Edgewise I believe
    J = Iy_local + Iz_local
    massmat[4,4] = rho_local*Iy_local
    massmat[5,5] = rho_local*Iz_local
    massmat[6,6] = rho_local*J

    stiffmat[1,1] = G*Area #Iz_local #flap shear stiffness #TODO: ???? I'm not sure if this is right. 
    stiffmat[2,2] = G*Area #Iy_local #edge shear stiffness
    stiffmat[3,3] = E*Area
    stiffmat[4,4] = E*Iy_local
    stiffmat[5,5] = E*Iz_local
    stiffmat[6,6] = G*J

    bdblade.nodes[i] = of.makenode(bdfracs[i], stiffmat, massmat)
end


####### Write files to disk
of.write_adfile(adfile, "NREL5MW_ADfile.dat";outputpath=path)

of.write_edfile(edfile, "NREL5MW_EDfile.dat"; outputpath=path) 
of.write_edblade(edblade, "NREL5MW_EDblade.dat"; outputpath=path)

of.write_bdblade(bdblade, "NREL5MW_BDblade.dat"; outputpath=path) #Todo: The solution won't converge with these values. :| I wonder if the blade is just to darn heavy. 





nothing