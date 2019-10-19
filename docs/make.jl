using SemiLagrangian
using Documenter

makedocs(modules=[SemiLagrangian],
         doctest = false,
         format = :html,
         sitename = "SemiLagrangian.jl",
         pages = ["Documentation" => "index.md",
                  "Functions"     => "functions.md",
                  "Contents"      => "contents.md"])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/JuliaVlasov/SemiLagrangian.jl.git",
 )
