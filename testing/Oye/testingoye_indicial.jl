using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, OpenFASTsr, LaTeXStrings

dsm = DynamicStallModels
of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)

c = 0.55

M = 0.11
a = 343.0
Vrel = M*a #60


dsmodel = Oye(Indicial(), 1, 2, 4.0)



af = dsm.make_airfoil(polar_0015, dsmodel, c; sfun=dsm.LSP())

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af




tvec = range(0, 1.0, 1000)
Uvec = Vrel.*ones(length(tvec))

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


states, loads = solve_indicial(airfoils, tvec, Uvec, alphavec)


stateplt = plot(tvec, states[:,1], leg=false, xaxis="time (s)", yaxis="f")



cn = loads[:,1]
if dsmodel.cflag == 2
    cn_static = af.cn.(alphavec)
else
    cn_static = af.cl.(alphavec)
end

cnplt = plot( xaxis="time (s)", yaxis=L"C_n", leg=:topright)
plot!(tvec, cn, lab="DSM")
plot!(tvec, cn_static, lab="Static")



cyclecnplt = plot(xaxis="Angle of Attack (deg)", yaxis=L"C_n", leg=:bottomright)
plot!(alphavec.*(180/pi), cn, lab="DSM")
plot!(alphavec.*(180/pi), cn_static, lab="Static")

display(cyclecnplt) 


nothing