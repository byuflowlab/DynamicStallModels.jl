using Test,DynamicStallModels,DelimitedFiles,OpenFASTTools
const repair_ds=DynamicStallModels
function repair_fixture(model=Oye(Discrete(),1,1,4.);sfun=(a->.25))
    a=[-.8,-.5,-.2,-.1,0.,.1,.2,.5,.8]
    cl=[-.4,-.7,-.6,-.4,0.,.4,.6,.7,.4];cd=fill(.02,9);cm=fill(.125,9)
    cn=cl.*cos.(a).+cd.*sin.(a);cc=cl.*sin.(a).-cd.*cos.(a)
    Airfoil(model,hcat(a,cl,cd,cm),Linear(a,cl),Linear(a,cd),Linear(a,cm),Linear(a,cn),Linear(a,cc),4.,4.,0.,[-.5,.5],[-pi/4,pi/4],pi/36,sfun)
end
@testset "Audited interface and numerical regressions" begin
    ds=repair_ds
    @testset "B01/B06 output-first state contract" begin
        af=repair_fixture();y=[2.,0.,.2,0.];p=[.5,.2];old=[.8];new=zeros(1)
        af(new,old,y,p,.1)
        @test old==[.8]
        @test af([1],y,p,.1)≈[.25+.75exp(-.1)]
        @test new≈[.25+.55exp(-.1)]
        @test af(old,y,p,.1)≈new
        @test old==[.8]
        afstep=repair_fixture(;sfun=a->a>.1 ? .25 : .8)
        s,l=solve_discrete([afstep],[0.,.1,.3],fill(2.,3),[0.,.2,.2];cvec=[.5])
        @test s[:,1]≈.25 .+.55exp.(-[0.,.1,.3])
        afc=repair_fixture(Oye(Continuous(),1,1,4.))
        result=[afc,afc]([1,2],[.8,.1],repeat(y,2),[.5,.2,1.,.2],0.)
        @test result≈[-.55,.075]
        @test size(l)==(3,3)
    end
    @testset "B02/B03/B04/B08/B10 polar constructors" begin
        polar=repair_fixture().polar;before=copy(polar)
        a=ds.make_airfoil(polar,ds.NoModel();sfun=ds.RSP())
        @test a.cm(.1)≈.125
        @test polar==before
        degrees=copy(polar);degrees[:,1].*=180/pi;degreecopy=copy(degrees)
        ad=ds.make_airfoil(degrees,ds.NoModel();sfun=ds.RSP(),radians=false)
        @test degrees==degreecopy
        @test ad.polar≈a.polar
        @test ad.alphasep≈a.alphasep
        @test ad.dcldalpha≈a.dcldalpha
        @test a.dcldalpha≈4.
        default=ds.make_airfoil(polar,ds.NoModel())
        @test all(isfinite,ds.separationpoint.(Ref(default),polar[:,1]))
        flat=copy(polar);flat[:,3].=0.
        defaultflat=ds.make_airfoil(flat,ds.NoModel())
        @test all(isfinite,ds.chordwiseseparationpoint.(Ref(defaultflat),flat[:,1]))
        simple=ds.make_simpleairfoil(polar[:,1:3],ds.NoModel())
        @test size(simple.polar,2)==4
        @test ds.update_airfoil(simple).cm(.1)==0.
        @test ds.make_simpleairfoil(polar,ds.NoModel()).cm(.1)≈.125
    end
    @testset "B07/B09/B13/B14/B15/B16" begin
        af=repair_fixture()
        bl=update_airfoil(af;sfun=ds.BLSP(fill(.05,4)))
        @test ds.separationpoint(bl,.5)≈.7
        @test ds.separationpoint(bl,-.5)≈.7
        @test_throws ArgumentError Linear([1.,0.],[1.,2.])
        @test_throws ArgumentError Linear([0.,0.],[1.,2.])
        @test_throws DimensionMismatch Linear([0.,1.,2.],[1.,2.])
        @test_throws ArgumentError Linear([0.,NaN],[1.,2.])
        @test_throws ArgumentError Linear([0.,1.],[1.,2.])(NaN)
        osp=update_airfoil(af;sfun=ds.OSP())
        @test isfinite(ds.separationpoint(osp,.1))
        @test ds.Clfs(0.,a->4a,4.,0.)==0.
        nm=update_airfoil(af;dsmodel=ds.NoModel())
        @test collect(ds.get_loads(nm.model,Float64[],[2.,0.,.1,0.],[.5,.2],nm))≈collect(ds.get_loads(nm.model,nm,Float64[],[2.,0.,.1,0.],[.5,.2]))
        @test ds.numberofparams(nm.model)==2
        @test ds.numberofloads(nm.model)==3
        @test_throws ArgumentError ds.NoModel(Continuous())
        @test all(s->isdefined(ds,s),names(ds))
    end
    @testset "B05/B20 ADG scalar/in-place/domain equivalence" begin
        repo=dirname(@__DIR__)
        af,xcp=OpenFASTTools.make_dsairfoil(OpenFASTTools.read_airfoilinput(joinpath(repo,"data","DU21_A17.dat")))
        y=[40.,0.,.1,0.];p=[1.,xcp];x,_=ds.initialize(af,[0.,.1],y,p)
        expected=copy(x);ds.update_states!(af.model,af,x,expected,y,p,.1)
        @test ds.update_states(af.model,af,x,y,p,.1)≈expected
        @test af(x,y,p,.1)≈expected
        @test ds.separationpoint(af,.1)==ds.separationpoint(af,.1,af.dcndalpha)
        @test x==ds.initialize(af,[0.,.1],y,p)[1]
        for speed in (-1.,0.,Inf,af.model.a)
            @test_throws ArgumentError solve_discrete([af],[0.,.1],fill(speed,2),[.1,.1])
        end
        for chord in (0.,-1.,Inf)
            @test_throws ArgumentError solve_discrete([af],[0.,.1],[40.,40.],[.1,.1];cvec=[chord])
        end
        @test_throws ArgumentError solve_discrete([af],[0.,.1],[40.,40.],[.1,NaN])
        @test_throws ArgumentError ds.update_states(af.model,af,x,y,p,0.)
        @test_throws DimensionMismatch ds.update_states(af.model,af,zeros(31),y,p,.1)
        @test_throws DimensionMismatch ds.BLADG_coefficients(af,zeros(31),y,p)
        for version in (1,2,4)
            m=BeddoesLeishman(Discrete(),version,af.model.A,af.model.b,af.model.T,af.model.Cn1,af.model.Cd0,af.model.Cm0,af.model.eta,af.model.zeta,af.model.a)
            @test_throws ArgumentError ds.initialize(update_airfoil(af;dsmodel=m),[0.,.1],y,p)
        end
    end
end
