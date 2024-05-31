using OpenFASTTools, FLOWMath, Plots

du21_a17 = of.read_airfoilinput("../../../data/DU21_A17.dat")

polar = hcat(du21_a17.aoa.*(pi/180), du21_a17.cl, du21_a17.cd, du21_a17.cm)
clfit = Akima(polar[:,1], polar[:,2])
cdfit = Akima(polar[:,1], polar[:,3])
cmfit = Akima(polar[:,1], polar[:,4])

alpha0 = du21_a17.alpha0*(pi/180)
aoa = polar[:,1]
cl = polar[:,2]
cdd = polar[:,3]
cd0 = cdfit(alpha0)

cn = @. cl*cos(aoa) + (cdd - cd0)*sin(aoa)
cc = @. cl*sin(aoa) - (cdd - cd0)*cos(aoa)


upplt = plot(xaxis="Angle of Attack (rads)", yaxis="Coefficient")
plot!(aoa, cl, lab="Lift")
plot!(aoa, cn, lab="Normal")
# display(upplt)

#=
- This looks like it makes sense. See once we get to high angles of attack, there should be a large separation bubble, which would mean a lot of drag. However, at this angle, the drag becomes the normal force... So it makes sense that it balloons back up after stall. 

- It appears that the lift curve slope is the same in the linear region. We don't really start to see differences until post stall. I thought about calculating the percent difference, but I think the difference will be less than the accuracy of approximating the lift curve slope. 
-> This tells me that I think it'd be safe to just use the lift curve slope for dcndalpha, and I think rather than computing the normal and tangential loads and storing them in the airfoil struct, I'll just compute them on the fly. Which is what Dr. Ning suggested, and OpenFAST does. It should be more efficient. 
=#


bplt = plot(xaxis="Angle of Attack (rads)", yaxis="Coefficient")
plot!(aoa, cdd, lab="Drag")
plot!(aoa, cc, lab="Tangent")
# display(bplt)s