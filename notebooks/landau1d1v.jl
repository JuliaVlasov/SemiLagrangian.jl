# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Julia 1.4.1
#     language: julia
#     name: julia-1.4
# ---

using LinearAlgebra
using ProgressMeter
using Plots

# +
import Base.Threads: @spawn, @sync, nthreads, threadid
abstract type InterpolationType end
abstract type AbstractAdvection end

include("../src/mesh.jl")
include("../src/geometry.jl")
include("../src/bspline_periodic.jl")
include("../src/advection.jl")

# +

eps    = 0.001
kx     = 0.5
degree = 3

xmin, xmax, nx =  0., 2π/kx, 32
vmin, vmax, nv = -6., 6., 64

mesh_x = UniformMesh( xmin, xmax, nx, endpoint = false )
mesh_v = UniformMesh( vmin, vmax, nv, endpoint = false )

advection_x! = Advection( mesh_x, Bspline(degree))
advection_v! = Advection( mesh_v, Bspline(degree))

ex  = zeros(Float64, nx)
rho = zeros(Float64, nx)

fxv = zeros(Float64, (nx, nv))
fvx = zeros(Float64, (nv, nx))

x = mesh_x.points
v = mesh_v.points

fxv .= (1 .+ eps .* cos.(kx .* x)) .* transpose(exp.(-.5*v.^2)) ./ sqrt(2π)

transpose!(fvx, fxv)

# +
tspan  = LinRange(0, 30, 600)
dt = 0.05

compute_charge!(rho, mesh_v, fvx)
compute_e!(ex, mesh_x, rho)

p = plot(layout=(1,2))
scatter!(p[1,1], x, rho, label=:computed, title="rho")
plot!(p[1,1], x, eps * cos.( kx .* x), label=:true)
scatter!(p[1,2], x , ex, label=:computed, title="Ex")
plot!(p[1,2], x , eps * sin.(kx .* x) / kx, label=:true)

# +
E = Float64[]

@showprogress 1 for t in tspan

    advection_x!(fxv, v, 0.5dt)
    transpose!(fvx, fxv)
    compute_charge!(rho, mesh_v, fvx)
    compute_e!(ex, mesh_x, rho)
    push!(E, 0.5 * log(sum(ex .* ex) * mesh_x.step))
    advection_v!(fvx, ex, dt)
    transpose!(fxv, fvx)
    advection_x!(fxv, v, 0.5dt)

end 

plot(tspan, E)
# -

