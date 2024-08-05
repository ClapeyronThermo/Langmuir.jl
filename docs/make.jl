using Documenter
using AdsorbedSolutionTheory


makedocs(;
  modules = AdsorbedSolutionTheory,
  authors = "Andres Riedemann, Vinicius Viena Santana and contributors",
  sitename = "AdsorbedSolutionTheory.jl",
  format = Documenter.HTML(;
    prettyurls = true,
    canonical = "https://ClapeyronThermo.github.io/AdsorbedSolutionTheory.jl/",
    assets = ["assets/themes/documenter-dark.css"],
  ),
  pages = [
    "Home" => "index.md",
    "Reference" => "reference.md"
  ],
)

deploydocs(; repo = "github.com/ClapeyronThermo/AdsorbedSolutionTheory.jl.git")