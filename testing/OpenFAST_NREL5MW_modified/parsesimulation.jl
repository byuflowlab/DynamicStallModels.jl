import OpenFASTsr
using Plots, FLOWMath, Statistics, DelimitedFiles

of = OpenFASTsr
fm = FLOWMath

path = dirname(@__FILE__)
cd(path)
# filepath = "../OpenFAST_files"

filename = "NREL5MW_input.out"

# outs = of.ReadOutput(filename, path)

fullout = readdlm(filename; skipstart=6)

names = fullout[1,:]
data = Float64.(fullout[3:end,:])

outs = Dict(names[i] => data[:,i] for i in 1:length(names))

adblade = of.read_adblade("NREL5MW_adblade.dat", path)

edblade = of.read_edblade("NREL5MW_edblade.dat", path)
edfile = of.read_edfile("NREL5MW_edfile.dat", path)

bdblade = of.read_bdblade("NREL5MW_BDblade.dat", path)

rvec_a = adblade.span .+ edfile.hubrad

rvec_elast = edblade.frac.*(edfile.tiprad-edfile.hubrad) .+ edfile.hubrad

nt = length(outs["Time"])
na = length(rvec_a)
ne = length(rvec_elast)
nb = length(bdblade.nodes)

rvec_beam = [bdblade.nodes[i].frac*(edfile.tiprad-edfile.hubrad) + edfile.hubrad for i = 1:nb]




# deflectionx_ed = zeros(nt, ne)
# deflectiony_ed = zeros(nt, ne)
# deflectionz_ed = zeros(nt, ne)

# for i = 1:ne      ####For use with ElastoDyn
#     if i<10
#         number = "00$i"
#     elseif i<100
#         number = "0$i"
#     else
#         number = "$i"
#     end
#     namex = "B1N"*number*"UxB"
#     namey = "B1N"*number*"UyB"
#     namez = "B1N"*number*"UzB"
#     deflectionx_ed[:,i] = outs[namex]
#     deflectiony_ed[:,i] = outs[namey]
#     deflectionz_ed[:,i] = outs[namez]
# end

deflectionx_beamdyn = zeros(nt, nb)
deflectiony_beamdyn = zeros(nt, nb)
deflectionz_beamdyn = zeros(nt, nb)

for i = 1:ne      ####For use with BeamDyn
    if i<10
        number = "00$i"
    elseif i<100
        number = "0$i"
    else
        number = "$i"
    end
    namex = "B1N"*number*"_TDxr"
    namey = "B1N"*number*"_TDyr"
    namez = "B1N"*number*"_TDzr"
    deflectionx_beamdyn[:,i] = outs[namex]
    deflectiony_beamdyn[:,i] = outs[namey]
    deflectionz_beamdyn[:,i] = outs[namez]
end


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


#Todo. I need to get the stiffnesses to match up between OpenFAST and GXBeam. I'm not sure what will be easier, using the stiffnesses calculated by Jonkman, or using my crap stiffnesses. 
# -> It would probably be best to use the stiffnesses from Jonkman, so that I can just use that as another validation case. 
# I suppose that I could take the coordinates file, and the chord length to get the shape, calculate Iyy and Izz. Then see how it compares to 
# -> I think it's going to be easiest to use the values that I came up with... because I know what the values should be. So I should be able to fill it in easy enough. 
# -> Okay, I think that I have them matching. I think at this point, if there is error, it is likely to be that I got Iyy and Izz mixed up. (That and in OpenFAST I'm simulating the cylinder region. )

tipdefplt = plot(xaxis="Time (s)", yaxis="Deflection (m)", legend=:topleft, title="Tip Deflection")
if @isdefined(deflectionx_beamdyn)
    plot!(tvec, deflectionx_beamdyn[:,end], lab="BeamDyn x")
    plot!(tvec, deflectiony_beamdyn[:,end], lab = "BeamDyn y")
    plot!(tvec, deflectionz_beamdyn[:,end], lab = "BeamDyn z")
end
# if @isdefined(deflectionx_ed)
#     plot!(tvec, deflectionx_ed[:,end], lab="ElastoDyn x", linestyle=:dash)
#     plot!(tvec, deflectiony_ed[:,end], lab="ElastoDyn y", linestyle=:dash)
#     plot!(tvec, deflectionz_ed[:,end], lab="ElastoDyn z", linestyle=:dash)
# end
display(tipdefplt)

#TODO: Regardless of if I'm off on Ixx and Iyy, it appears that my solution (using the fixed point solution) is off by an order of magnitude. ... I should check to see if the aerodynamic solution is off by some amount... 

forceplt = plot(xaxis="Radius (m)", yaxis="Distributed Load (N/m)", legend=:topleft)
plot!(rvec_a, fxmat[end-6000,:], lab="Fx")
plot!(rvec_a, fymat[end-6000,:], lab="Fy")
display(forceplt)

# forceanim = @animate for i in 1:nt
#     plot(xaxis="Radius (m)", yaxis="Distributed Load (N/m)", legend=:topleft, xlims = (0,65), ylims=(-300, 8500))
#     plot!(rvec_a, fxmat[i,:], lab="Fx")
#     plot!(rvec_a, fymat[i,:], lab="Fy")
# end
# gif(forceanim, "/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/openfast_simulation/figures/forceanim_fps1000.gif", fps = 1000)

# displacementanim = @animate for i in 1:nt
#     plot(xaxis="Radius (m)", yaxis="Displacement (m)", legend=:topleft, xlims=(0,65), ylims=(-1,5))
#     plot!(rvec_beam, deflectionx_beamdyn[i,:], lab="Ux")
#     plot!(rvec_beam, deflectiony_beamdyn[i,:], lab="Uy")
#     plot!(rvec_beam, deflectionz_beamdyn[i,:], lab="Uz")
# end
# gif(displacementanim, "/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/openfast_simulation/figures/displacementanim_fps1000.gif", fps = 1000)


#=
The gifs are a slow... so heads up. slow to run. Slow to display. slow to compile. 

The solution didn't change a whole ton when I converted from the stiffness and mass matrices from precomp to the ones I prepared.... which is.... odd. I would have expected it to change. I guess that just suggests that my model isn't incredibly far from reality? That or I have something wrong. So that's comforting. 

=#