using Documenter
using Langmuir


makedocs(;
  modules = Langmuir,
  checkdocs=:none,
  authors = "Andres Riedemann, Vinicius Viena Santana and contributors",
  sitename = "Langmuir.jl",
  format = Documenter.HTML(;
    prettyurls = true,
    canonical = "https://ClapeyronThermo.github.io/Langmuir.jl/",
    assets = ["assets/logo.ico"],
  ),
  pages = [
    "Home" => "index.md",
    "Background" => "tutorials/background.md",
    "Getting Started" => "tutorials/getting_started.md",
    "Tutorials" => Any["tutorials/tutorial.md"],
    "Supported models" => "models/models.md",
    "Reference" => "reference.md",
  ],
)

deploydocs(; repo = "github.com/ClapeyronThermo/Langmuir.jl.git", push_preview = true)