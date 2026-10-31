# Separate from unit tests: includes plotting and an ODE solver. Run with
# julia --project=examples test/run_documentation_examples.jl
using DynamicStallModels, Test
repo=dirname(@__DIR__)
@assert realpath(pathof(DynamicStallModels))==realpath(joinpath(repo,"src","DynamicStallModels.jl"))
cd(repo)
@testset "Runnable documentation" begin
    for page in ("gettingstarted.md","beddoesleishman_example.md","Oye_example.md")
        @testset "$page" begin
            scope=Module(gensym(:DocumentationExample))
            source=read(joinpath(repo,"docs","src",page),String)
            blocks=collect(eachmatch(r"(?ms)^```julia\s*\n(.*?)^```",source))
            @test !isempty(blocks)
            for (i,block) in enumerate(blocks)
                Base.include_string(scope,block.captures[1],page*":block"*string(i))
            end
            if haskey(ENV,"DOCUMENTATION_PLOT_DIR") && isdefined(scope,:plot_path)
                mkpath(ENV["DOCUMENTATION_PLOT_DIR"])
                cp(getfield(scope,:plot_path),joinpath(ENV["DOCUMENTATION_PLOT_DIR"],page*".png");force=true)
            end
            if page=="Oye_example.md"
                @test maximum(abs.(getfield(scope,:states)[:,1].-vec(Array(getfield(scope,:solution)))))<0.02
                @test all(f -> -1e-8<=f<=1+1e-8, Array(getfield(scope,:solution)))
            else
                @test size(getfield(scope,:loads))==(100,3)
                @test getfield(scope,:states)[:,1]==getfield(scope,:alpha)
            end
        end
    end
end
