using Plots

"""
11/5/21 Adam Cardoza
Testing the myrk4 method that I wrote to have a "stupid" ODE solver on hand. 
"""

include("testingdynamicstallmethods.jl")

function states!(dx, x, p, t)
    dx[1] = x[2]
    dx[2] = -4*x[1]
end

function solution(t)
    return cos(2*t)
end

x0 = [1.0, 0.0]
p = [1]
rng = (0.0, 10.0)

tvec, xvec = myrk4(states!, x0, rng, p; numsteps=200)

treal = collect(range(rng[1], rng[2], length=200))
xreal = solution.(treal)

plt = plot(legend=:topright)
plot!(treal, xreal, lab="Real")
scatter!(tvec, xvec[:,1], lab="RK4")
display(plt)

### My conclusion is that this is a good simple method. for a range of 10, the solution converged with about 25 steps, which in my opinion, is pretty dang good. It might not be the fastest solution method... and who knows what it is weak to... but it should work. 