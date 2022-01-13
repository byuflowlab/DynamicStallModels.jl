using Plots, DelimitedFiles, Xfoil, FLOWMath

polar = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Leishman1989/leishman1989staticdata_fig1.csv", ',')

# extract geometry
x = Float64[]
y = Float64[]

f = open("/Users/adamcardoza/Box/research/FLOW/bladeopt/data/coordinates/naca0012.dat", "r")

for line in eachline(f)
    entries = split(chomp(line))
    push!(x, parse(Float64, entries[1]))
    push!(y, parse(Float64, entries[2]))
end

close(f)

# set operating conditions
c = 0.1
M = 0.379
a = 343.0
v = M*a
rho = 1.225
mu = 18.03e-6
# re = rho*v*c/mu
alpha = -15:1:30
# re = 3.83e6
re = 4e6

c_l, c_d, c_dp, c_m, converged = Xfoil.alpha_sweep(x, y, alpha, re, mach=0.383, iter=200, zeroinit=false, printdata=false)


liftplt = plot(polar[:,1], polar[:,2], leg=:topleft, lab="Polar")
plot!(alpha, c_l, lab="Xfoil")
xlabel!("Angle of Attack (degrees)")
ylabel!("Lift Coefficient")
display(liftplt)

dragplt = plot(leg=:topleft)
plot!(alpha, c_d, lab="Cd")
plot!(alpha, c_dp, lab="Cdp")
xlabel!("Angle of Attack (degrees)")
ylabel!("Drag Coeffficient")
display(dragplt)

nothing