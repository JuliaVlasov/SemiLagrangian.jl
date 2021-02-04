push!(LOAD_PATH,"../src/")

using Documenter
using SemiLagrangian
using Plots 

ENV["GKSwstype"] = "100"

makedocs(modules=[SemiLagrangian],
         sitename = "SemiLagrangian.jl",
         authors = "Yves Mocquard, Pierre Navaro and Nicolas Crouseilles",
         format=Documenter.HTML(;
         prettyurls=get(ENV, "CI", "false") == "true",
         canonical="https://juliavlasov.github.io/SemiLagrangian.jl",
         assets=String[],
         ),
         doctest = false,
         pages = ["Home" => "index.md",
                  "Quickstart" => "quickstart.md",
                  "Functions"     => "functions.md",
                  "Contents"      => "contents.md"])

deploydocs(;
    branch = "gh-pages",
    devbranch = "master",
    repo   = "github.com/JuliaVlasov/SemiLagrangian.jl"
)
