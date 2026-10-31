using Test,DynamicStallModels
module RepairedLegacyBL
include("../src/BeddoesLeishman/indicial/BeddoesLeishman.jl")
end
module RepairedLegacyLarsen
using DynamicStallModels: DSModel,DEType,Airfoil
include("../src/Larsen.jl")
end
struct RepairSolution
    t::Vector{Float64}
    values::Matrix{Float64}
end
Base.Array(sol::RepairSolution)=sol.values
@testset "Explicit model parameters and legacy numerical repairs" begin
    ds=DynamicStallModels
    @testset "B17/B18 homogeneous relaxation and reduced distance" begin
        legacy=RepairedLegacyBL
        for step in (0.,.01,1.,10.),time in (.1,1.,10.)
            expected=exp(-step/time)
            @test legacy.CvNfun(1.,.5,.5,step,time)≈expected
            @test legacy.Dffun(1.,step,time,.5,.5)≈expected
            @test legacy.DPfun(1.,step,time,.5,.5)≈expected
        end
        @test legacy.delSfun(1.,20.,40.,1.1,1.)≈6.
        @test legacy.delSfun(2.,40.,40.,2.,1.)≈40.
        @test legacy.delSfun(1.,0.,0.,2.,1.)==0.
        @test_throws ArgumentError legacy.delSfun(0.,20.,40.,1.1,1.)
        @test_throws ArgumentError legacy.CvNfun(1.,0.,0.,1.,0.)
    end
    @testset "B19 Larsen explicit closure inputs" begin
        m=RepairedLegacyLarsen
        p=(CL0=t->1.,CL0dot=t->2.,alpha=t->.1,alphadot=t->.01,alphav=.2,w=ones(4),A=[.3,.7],theta=t->.4)
        du=zeros(4);x=zeros(4)
        m.states!(du,x,p,0.)
        @test du≈[.6,1.4,.4,0.]
        @test x==zeros(4)
        @test_throws ArgumentError m.states!(du,x,collect(values(p))[1:6],0.)
    end
    @testset "Risø wrappers, B22/B23 residuals and station extraction" begin
        af=repair_fixture(;sfun=ds.RSP())
        A=[.3,.7];b=[.14,.53];T=[1.7,3.]
        model=Riso([af,af];A,b,T)
        @test_throws ArgumentError Riso([af])
        @test_throws ArgumentError Riso([af];A,b,T,detype=Discrete())
        p=[2.,0.,.1,0.,.5,1.]
        one=[.03,.07,.4,ds.separationpoint(af,.1)]
        x=vcat(one,one);dx=zeros(8)
        model(dx,x,p,0.)
        @test maximum(abs,dx)<1e-14
        @test model(x,p,0.)≈dx
        @test x==vcat(one,one)
        changed=vcat(one,[.01,.02,.2,.3]);rates=model(changed,p,0.)
        @test rates[1:4]!=rates[5:8]
        @test ds.riso_residual(zeros(4),one,p[1:4],.5,0.,af;A,b,T)≈zeros(4) atol=1e-14
        @test_throws ArgumentError ds.riso_residual(zeros(4),one,p[1:4],.5,0.,af)
        @test_throws ArgumentError ds.riso_residual(zeros(4),one,[2.,0.,1.,0.,.1,0.],.5,0.,af;A,b,T)
        explicit=[2.,0.,.1,0.,.5,A[1],A[2],b[1],b[2],af.dcldalpha,af.alpha0,T[1],T[2],af.cl,af.alphasep[1],af.alphasep[2]]
        res=zeros(4);derivative=[1.,2.,3.,4.]
        ds.static_riso_dae!(res,derivative,one,explicit,0.,af)
        @test res≈derivative atol=1e-14
        @test derivative==[1.,2.,3.,4.]
        @test_throws ArgumentError ds.static_riso_dae!(res,derivative,one,explicit,0.)
        sol=RepairSolution([0.,.1],hcat(x,changed))
        cl,cd,t=parsesolution(model,sol,p)
        @test size(cl)==size(cd)==(2,2)
        @test cl[2,1]!=cl[2,2]
        @test cl[2,2]≈ds.riso_coefficients(changed[5:8],2.,.1,0.,1.,af;A)[1]
        @test all(isfinite,cl) && all(isfinite,cd)
    end
end
