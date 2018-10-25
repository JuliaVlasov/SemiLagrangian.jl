push!(LOAD_PATH,"../src/")

using SplittingOperators
using Documenter
using Plots # to not capture precompilation output

makedocs(modules=[SemiLagrangian],
         doctest = false,
         format = :html,
         sitename = "SemiLagrangian.jl",
         pages = [ "Semi-Lagrangian" => "bsl.md",
		   "Advection functions" => "advections.md",
		   "Contents" => "contents.md"])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/JuliaVlasov/SemiLagrangian.jl.git",
    julia  = "1.0",
    osname = "osx"
 )
