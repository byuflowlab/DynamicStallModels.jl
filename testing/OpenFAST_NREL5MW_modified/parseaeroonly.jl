import OpenFASTsr
using Plots, FLOWMath, Statistics, DelimitedFiles

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

edblade = of.read_edblade("NREL5MW_edblade.dat", path)
edfile = of.read_edfile("NREL5MW_edfile.dat", path)


rvec_a = adblade.span .+ edfile.hubrad

nt = length(outs["Time"])
na = length(rvec_a)




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



forceplt = plot(xaxis="Radius (m)", yaxis="Distributed Load (N/m)", legend=:topleft)
plot!(rvec_a, fxmat[end,:], lab="Fx")
plot!(rvec_a, fymat[end,:], lab="Fy")
display(forceplt)


tipplt = plot(xaxis="Time (s)", yaxis="Distributed Load (N/m)", title="Tip Loads")
plot!(tvec, fxmat[:,end-1], lab="Fx")
plot!(tvec, fymat[:,end-1], lab="Fy")
display(tipplt)

