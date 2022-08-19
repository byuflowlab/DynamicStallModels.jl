using Documenter, DynamicStallModels

makedocs(
        modules = [DynamicStallModels],
        format = Documenter.HTML(),
        pages = [
    "Introduction" => "index.md",
    "Getting Started" => "gettingstarted.md",
    "Examples" => ["riso_example.md"],
    "Developers" => ["codeoutline.md", "riso_theory.md"],
    "API Reference" => "apireference.md"],
        repo="https://github.com/byuflowlab/DynamicStallModels.jl/blob/{commit}{path}#L{line}",
        sitename="DynamicStallModels.jl",
        authors = "Adam Cardoza <adam.cardoza.online@gmail.com>",
        )

deploydocs(
    repo = "github.com/byuflowlab/DynamicStallModels.jl.git",
)