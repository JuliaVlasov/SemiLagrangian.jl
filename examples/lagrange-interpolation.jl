# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Julia 1.8.5
#     language: julia
#     name: julia-1.8
# ---

# # Lagrange Interpolation

using SemiLagrangian
using Plots
using Test

order = 5
lag = Lagrange(order, Float64)

f(x) = (x+2) * (3x+2) * (2x-1) * (x - 2) * (3x+1)

xi = LinRange(-2, 2, 20)
fi = Float64[]
for x in xi
    res = sum(lag.tabfct[j+1](x) * f(j - order ÷ 2) for j = 0:order)
    push!(fi, res)
end

plot( f, LinRange(-2.1, 2.1, 100))
scatter!(xi, fi)

g(x, y) = exp(-(x^2+y^2))

# +
order = 3
lag = Lagrange(order, Float64)
nx, ny = 100, 100
xi = LinRange(-4, 4, nx)
yi = LinRange(-4, 4, ny)
gi = zeros(nx, ny)

for i in eachindex(xi), j in eachindex(yi)
    dec = order ÷ 2
    res = 0.0
    for l = 0:order, k = 0:order
        res += lag.tabfct[l+1](xi[i]) * lag.tabfct[k+1](yi[j]) * g(l - dec, k - dec)
    end
    gi[i, j] = res
    
end

# -

contour(xi, yi, g)

contour(xi, yi, gi)


