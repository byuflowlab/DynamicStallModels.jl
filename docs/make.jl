using Documenter, DynamicStallModels

makedocs(sitename="DynamicStallModels.jl", pages = [
    "Introduction" => "index.md",
    "Getting Started" => "gettingstarted.md",
    "Examples" => ["riso_example.md"],
    "Developers" => ["codeoutline.md", "riso_theory.md"],
    "API Reference" => "apireference.md"
])

deploydocs(
    repo = "github.com/byuflowlab/DynamicStallModels.jl.git",
)