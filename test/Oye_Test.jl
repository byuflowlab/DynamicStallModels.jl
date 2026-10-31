# Numerical verification, not an experimental validation of dynamic stall.
using DynamicStallModels, Test

function oye_test_airfoil(model; sfun=(alpha -> 0.25), lift=nothing)
    a = [-0.7,-0.5,-0.2,-0.1,0.0,0.1,0.2,0.5,0.7]
    cl = isnothing(lift) ? [-0.4,-0.7,-0.6,-0.4,0.0,0.4,0.6,0.7,0.4] : lift
    cd = fill(0.02,length(a)); cm = zeros(length(a))
    cn = cl.*cos.(a).+cd.*sin.(a); cc = cl.*sin.(a).-cd.*cos.(a)
    return Airfoil(model,hcat(a,cl,cd,cm),Linear(a,cl),Linear(a,cd),Linear(a,cm),
        Linear(a,cn),Linear(a,cc),4.0,4.0,0.0,[-0.5,0.5],[-pi/4,pi/4],pi/36,sfun)
end

struct OyeTestSolution
    t::Vector{Float64}
    states::Matrix{Float64}
end
Base.Array(sol::OyeTestSolution) = sol.states

@testset "Øye" begin
    ds = DynamicStallModels
    @testset "Supported configuration and domain" begin
        for flag in (-1,0,2,3)
            @test_throws ArgumentError Oye(Discrete(),flag,1,4.0)
        end
        for version in (0,3)
            @test_throws ArgumentError Oye(Discrete(),1,version,4.0)
        end
        for A in (-1.0,0.0,NaN,Inf)
            @test_throws ArgumentError Oye(Discrete(),1,1,A)
        end
        af=oye_test_airfoil(Oye(Discrete(),1,1,4.0))
        @test ds.numberofstates(af.model)==1
        @test ds.numberofloads(af.model)==3
        @test ds.update_states(af.model,af,[0.8],[0.0,0.0,0.2,0.0],[0.5,0.2],1.0)==[0.8]
        @test ds.update_states(af.model,af,[0.8],[2.0,0.0,0.2,0.0],[0.5,0.2],0.0)==[0.8]
        @test_throws ArgumentError ds.update_states(af.model,af,[0.8],[-1.0,0.0,0.2,0.0],[0.5,0.2],0.1)
        @test_throws ArgumentError ds.update_states(af.model,af,[0.8],[2.0,0.0,0.2,0.0],[0.0,0.2],0.1)
        @test_throws ArgumentError ds.update_states(af.model,af,[0.8],[2.0,0.0,0.2,0.0],[0.5,0.2],-0.1)
    end
    @testset "Closed-form separation step and solver contract" begin
        af=oye_test_airfoil(Oye(Discrete(),1,1,4.0);sfun=(a -> a>0.1 ? 0.25 : 0.8))
        t=[0.0,0.1,0.2,0.3];U=fill(2.0,4);alpha=[0.0,0.2,0.2,0.2]
        states,loads=solve_discrete([af],t,U,alpha;cvec=[0.5])
        @test states[:,1] ≈ 0.25 .+ 0.55.*exp.(-t) # tau=A*c/U=1 second
        @test size(states)==(4,1) && size(loads)==(4,3)
        @test all(isfinite,loads) && all(x -> 0<=x<=1,states)
        @test all(iszero,loads[:,2:3]) # lift-only model contract
        @test solve_indicial === solve_discrete
        @test solve_indicial([af],t,U,alpha;cvec=[0.5])==(states,loads)
        @test solve_discrete(view([af],:),t,U,alpha;cvec=[0.5])==(states,loads)
        one_state,one_load=solve_discrete([af],[0.0],[2.0],[0.0])
        @test size(one_state)==(1,1) && size(one_load)==(1,3)
        @test_throws ArgumentError solve_discrete([af],[0.0,0.0],[2.0,2.0],[0.0,0.0])
        @test_throws DimensionMismatch solve_discrete([af],t,U[1:2],alpha)
        @test_throws DimensionMismatch solve_discrete([af],t,U,alpha;cvec=Float64[])
        afc=oye_test_airfoil(Oye(Continuous(),1,1,4.0))
        @test_throws ArgumentError solve_discrete([afc],t,U,alpha)
    end
    @testset "Finite attached limit and static polar recovery" begin
        for version in (1,2)
            af=oye_test_airfoil(Oye(Discrete(),1,version,4.0);
                sfun=version==1 ? ds.RSP() : ds.LSP())
            for a in (-0.6,-0.2,0.0,0.1,0.2,0.5,0.6)
                states,loads=ds.initialize(af,[0.0],[10.0,0.0,a,0.0],[0.55,0.2])
                @test all(isfinite,states) && all(isfinite,loads)
                @test 0<=states[1]<=1
                # Exact equilibrium recovery where the selected separation
                # convention is consistent; LSP intentionally has a .35 rule.
                if version==1
                    @test loads[1] ≈ af.cl(a) atol=1e-12
                end
            end
        end
        attached=oye_test_airfoil(Oye(Discrete(),1,1,4.0);sfun=(a->1.0))
        @test ds.cl_fullysep_hansen(attached,0.0)==0.0
        @test ds.cl_fullysep_hansen(attached,0.1) ≈ 0.2
        # A low static coefficient can give negative raw LSP. It must be a
        # numeric zero, not `nothing` from the historical bare return.
        low=oye_test_airfoil(Oye(Discrete(),1,2,4.0);sfun=ds.LSP(),
            lift=[-0.1,-0.1,-0.01,-0.005,0.0,0.005,0.01,0.1,0.1])
        @test ds.separationpoint(low,0.2)==0.0
        @test ds.oye_attachment(oye_test_airfoil(low.model;sfun=(a->-0.1)),0.2)==0.0
        @test ds.oye_attachment(oye_test_airfoil(low.model;sfun=(a->1.1)),0.2)==1.0
    end
    @testset "Hermite endpoint values and slopes" begin
        af=oye_test_airfoil(Oye(Discrete(),1,2,4.0);sfun=ds.LSP())
        h=1e-6;f(a)=ds.cl_fullysep_faber(af,a)
        @test f(0.0)==af.cl(0.0)
        @test f(0.5)==af.cl(0.5)
        @test (f(h)-f(0.0))/h ≈ af.dcldalpha/2 rtol=1e-4
        @test (f(0.5)-f(0.5-h))/h ≈ af.dcldalpha/12 rtol=1e-4
    end
    @testset "Continuous dispatch, geometry and parsing" begin
        a1=oye_test_airfoil(Oye(Continuous(),1,1,4.0))
        a2=oye_test_airfoil(Oye(Continuous(),1,1,4.0);sfun=(a->0.75))
        y=[t->2.0,0.0,t->0.2,0.0,t->2.0,0.0,t->0.2,0.0]
        x=[0.8,0.1];dx=zeros(2)
        prob=ODEProblem([a1,a2],x,(0.0,1.0),y;cvec=[0.5,1.0])
        prob.f(dx,x,prob.p,0.2)
        @test dx ≈ [-0.55,0.325]
        @test x==[0.8,0.1] # rates must not overwrite the state
        @test_throws DimensionMismatch ODEProblem([a1],x,(0.0,1.0),y[1:4])
        @test_throws DimensionMismatch ODEProblem([a1],[0.8],(0.0,1.0),y)
        fake=OyeTestSolution([0.0,0.5],reshape([0.8,0.25],1,2))
        parsed=parsesolution(a1.model,a1,fake,y[1:4],[0.5,0.2])
        @test parsed[1,:]==[0.2,0.2]
        @test parsed[2,end] ≈ a1.cl(0.2)
        @test parsed[2,1] ≈ 0.8*0.8+0.2*((0.6-0.25*0.8)/0.75)
    end
end
