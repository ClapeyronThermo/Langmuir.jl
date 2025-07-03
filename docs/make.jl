using Documenter
using Langmuir
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "refs.bib"))

makedocs(;
  modules = Langmuir,
  checkdocs=:none,
  authors = "Andres Riedemann, Vinicius Viena Santana and contributors",
  sitename = "Langmuir.jl",
  plugins = [bib],
  format = Documenter.HTML(;
    prettyurls = true,
    canonical = "https://ClapeyronThermo.github.io/Langmuir.jl/",
    assets = ["assets/logo.ico"],
  ),
  pages = [
    "Home" => "index.md",
    "Background" => "tutorials/background.md",
    "Getting Started" => "tutorials/getting_started.md",
    "Tutorials" => Any["tutorials/tutorial.md", "tutorials/isosteric_heat.md", "tutorials/multicomponent.md"],
    "Supported models" => "models/models.md",
    "Reference" => "reference.md",
    "Bibliography" => "references.md",
  ],
)

deploydocs(; repo = "github.com/ClapeyronThermo/Langmuir.jl.git", push_preview = true)