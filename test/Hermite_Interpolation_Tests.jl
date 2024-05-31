using DynamicStallModels, DelimitedFiles, Plots, FLOWMath, OpenFASTTools, LaTeXStrings
using Test

dsm = DynamicStallModels
of = OpenFASTTools

path = dirname(@__FILE__)
cd(path)

c = 0.55

M = 0.0791
a = 343.0
Vrel = M*a #60

# polar = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/learning/exploring/exampledata/NACA4412.dat", '\t'; skipstart=3) 
# af = airfoil(polar) #Todo: This constructor is broken.

#du21_a17 = of.read_airfoilinput("../../data/airfoils/DU40_A17.dat") 
#af = of.make_dsairfoil(du21_a17) #Todo: I think this polar might be my problem. I should try a different polar.... which means that I need to fix the constructor. :| 

# af = update_airfoil(af, A=[4.0], dcndalpha=6.320368333107256, alpha0=-0.0033903071711640564)

af = airfoil(polar_0015; A=[12.3324])

af_new = update_airfoil(af)

airfoils = Array{Airfoil, 1}(undef, 1)
airfoils[1] = af_new


dsmodel = Oye(Indicial(), 1, airfoils, 1, 2)


tvec = range(0, 1.0, 1000) #0:0.001:0.05
Uvec = Vrel.*ones(length(tvec))

function alpha(t)
    c = 0.55
    M = 0.0791
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

states, loads = solve_indicial(dsmodel, [c], tvec, Uvec, alphavec)


"""
stateplt = plot(tvec, states[:,1], leg=false, xaxis="time (s)", yaxis="f")
# display(stateplt)


cn = loads[:,1]
cn_static = af.cn.(alphavec)

cnplt = plot( xaxis="time (s)", yaxis=L"C_l", leg=:topright)
plot!(tvec, cn, lab="DSM")
plot!(tvec, cn_static, lab="Static")
# display(cnplt) 


cyclecnplt = plot(xaxis="Angle of Attack (deg)", yaxis=L"C_l", leg=:topright)
plot!((alphavec.*(180/pi)), cn, lab="DSM")
plot!(alphavec.*(180/pi), cn_static, lab="Static")
display(cyclecnplt) 

#=
I think that I have the models matching, the difference is that I'm using Cn in this model and Cl in the other. 
=#
"""

#additional change make the equality statement at the beginning an equal to as well
@testset "Hermite Testing" begin
    function cl_fullysep_faber(cl, cl_sep, dcldalpha, alpha, alpha0, alpha_sep)
        if alpha0 < alpha < alpha_sep
            #Hermite Interpolation as in Faber 2018
            #TODO: Needs to be extended to work for angles less than alpha0, it should have a reflection to what is done here. 
            t0 = (alpha - alpha0)/(alpha_sep - alpha0) #Faber 2018 EQ A.1a
            #As alpha -> alpha_sep this will diverge. So... what keeps it from diverging? 
            t1 = (alpha - alpha_sep)/(alpha_sep - alpha0) #Faber 2018 EQ A.1b
            term1 = (alpha_sep - alpha0)*dcldalpha*(1 + t0*(7*t1/6 - 1))/2
            term2 = cl_sep*t0*(1-( 2*t1))
            # return t0*(term1+term2) #Faber 2018 EQ A.4
            clfs = t0*(term1+term2) #Faber 2018 EQ A.4
    
            # if clfs >10
            #     @show t0, t1, term1, term2
            # end
    
            if t0>10
                @show alpha, alpha0, alpha_sep
            end
    
            return clfs
        else
            return cl
        end
    end

    polar = polar_0015
    cl = Akima(polar[:,1], polar[:,2])
    _, minclidx = findmin(polar[:,2])
    _, maxclidx = findmax(polar[:,2])
    alphasep = [-32*pi/180, 32*pi/180]
    alpha0, _ = brent(cl, -0.25, 0.25)

    #@testset "clfs" for i in 1:length(polar[:,1])
        #@test cl_fullysep_faber(cl, cl(alphasep[2]), 2*pi, polar[i,1], alpha0, alphasep[2]) >= 0
    #end
    @testset "Hermite Test" for i in 1:length(polar[:,1])
       @test cl_fullysep_faber(cl, cl(alphasep[2]), 2*pi, polar[i,1], alpha0, alphasep[2]) < cl(polar[i,1])
    end

end



nothing
