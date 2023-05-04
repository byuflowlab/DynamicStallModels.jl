using DynamicStallModels, DelimitedFiles
dsm = DynamicStallModels
polar_0015 = readdlm("NACA_0015_Faber.csv" , ',')
dsmodel = Oye(Indicial(), 1, 2, 4.0)

c = 0.55
M = 0.11
a = 343.0
Vrel = M*a 

af = dsm.make_airfoil(polar_0015, dsmodel, c; sfun=dsm.LSP())
airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af


tvec = range(0, 1.0, 1000)
function alpha(t)
    c = 0.55
    M = 0.11
    a = 343.0
    shift = 10.0
    amp = 10.0
    k = 0.051

    v = M*a
    omega = k*2*v/c

    alf = shift + amp*sin(omega*t)
    return alf*(pi/180)
end

alphavec = alpha.(tvec)
Uvec = Vrel.*ones(length(tvec))


states, loads = solve_indicial(airfoils, tvec, Uvec, alphavec)


cn = loads[:,1]
if dsmodel.cflag == 2
    cn_static = af.cn.(alphavec)
else
    cn_static = af.cl.(alphavec)
end


using Plots, LaTeXStrings
stateplt = plot(tvec, states[:,1], leg=false, xaxis="time (s)", yaxis="f")
cyclecnplt = plot(xaxis="Angle of Attack (deg)", yaxis=L"C_L", leg=:bottomright)
plot!(alphavec.*(180/pi), cn, lab="DSM")
plot!(alphavec.*(180/pi), cn_static, lab="Static")
display(cyclecnplt) 