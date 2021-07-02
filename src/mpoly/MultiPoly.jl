#module MultiPoly

using DataStructures
using Combinatorics
using LinearAlgebra

export
    MPoly, terms, vars, nvars, generators, generator, monomials, deg,
    PolyUnion, newexps, oldexps,
    harmonicpolynomialdim, laplaceharmonicpol,
    evaluate, evaluate_basis, diff, integrate

import Base: zero, one,
    show, showcompact, print, length, endof, getindex, setindex!, copy, promote_rule, convert, eltype,
    *, /, //, -, +, ==, ^, divrem, conj, rem, real, imag, diff

include("mpoly.jl")
include("mpolyarithmetic.jl")
include("mpolyprinting.jl")
include("polyunion.jl")
include("mpolyevaluation.jl")
include("mpolycalculus.jl")

#end # module
