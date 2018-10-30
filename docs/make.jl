push!(LOAD_PATH,"../src/")

using SemiLagrangian
using Documenter

makedocs(modules=[SemiLagrangian],
         doctest = false,
         format = :html,
         sitename = "SemiLagrangian.jl",
         pages = ["Documentation" => "index.md",
                  "Advections"    => "advections.md",
                  "BSL"           => "bsl.md",
                  "Contents"      => "contents.md"]

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/JuliaVlasov/SemiLagrangian.jl.git",
 )
