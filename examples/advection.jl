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

# # Advection

using LinearAlgebra
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

g(x, y) = f(x) * f(y)

# +
order = 5
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

@test f.(xi) .* f.(yi') ≈ gi

fct1(x) = cos(2π * x + 0.25)
fct2(x) = exp(-(cos(2π * x + 0.25) - 1)^2)
fct3(x) = exp(-((x - 0.4)/0.1)^2)

nx = 128
mesh = UniformMesh(0., 1., nx)
fp = fct3.(mesh.points)
fi = copy(fp)
order = 5
interp = Lagrange(order, Float64)


# +
function advection!(fp, fi, mesh, interp, v, dt)

    dec = - v * dt / mesh.step
    decint = floor(Int, dec)
    value = dec - decint
    
    if get_order(interp) % 2 == 0 && value > 0.5
        value -= 1
        decint += 1
    end
        
    precal =  [fct(value) for fct in interp.tabfct]    
    
    interpolate!(fp, fi, decint, precal, interp)
    
end
# -

v, dt = 1.0, 1.0
nsteps = 100
dt = 0.2 / nsteps
order = 5
fp = fct3.(mesh.points)
fi = copy(fp)
plot(mesh.points, fp, xlims=(0,1))
interp = Lagrange(order, Float64)
interp = BSplineLU(order, nx, Float64)
interp = BSplineFFT(order, nx, Float64)
for i in 1:nsteps
    advection!(fp, fi, mesh, interp, v, dt)
    fi .= fp
end
plot!(mesh.points, fp, xlims=(0,1))


