# SemiLagrangian.jl


![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.org/JuliaVlasov/SemiLagrangian.jl.svg?branch=master)](https://travis-ci.org/JuliaVlasov/SemiLagrangian.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliavlasov.github.io/SemiLagrangian.jl/dev)
[![codecov](https://codecov.io/gh/JuliaVlasov/SemiLagrangian.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaVlasov/SemiLagrangian.jl)

Advection implementation using Semi-Lagrangian numerical method

**NOTE: This package is still very much under development and is not fully tested.**

## Installation

```julia-repl
julia>]
pkg> registry add https://github.com/juliavlasov/Registry
pkg> add SemiLagrangian
julia> using SemiLagrangian
```

Some examples are available in notebooks directory.

## Fortran code

The original Fortran code was used for the paper

"Charge conserving grid based methods for the Vlasov-Maxwell equations, Comptes Rendus de Mécanique (Theoritical and numerical approaches for Vlasov-Maxwell equations), 342, pp. 636-646 (2014)."

by N. Crouseilles, P. Navaro, E. Sonnendrücker.

http://people.rennes.inria.fr/Nicolas.Crouseilles/charge-cons-eulerian.pdf

https://doi.org/10.1016/j.crme.2014.06.012
