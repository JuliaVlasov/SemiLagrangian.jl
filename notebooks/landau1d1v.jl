# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.2.0
#     language: julia
#     name: julia-1.2
# ---

# +

abstract type InterpolationType end
abstract type AbstractAdvection end

include("../src/mesh.jl")
include("../src/geometry.jl")
include("../src/bspline_periodic.jl")
include("../src/splinenn.jl")
include("../src/splinepp.jl")
include("../src/banded_matrix.jl")
include("../src/spline_1d.jl")
include("../src/spline_interpolator_1d.jl")
include("../src/advection.jl")



# +
using LinearAlgebra, Plots, ProgressMeter

eps    = 0.001
kx     = 0.5
degree = 5

xmin, xmax, nx =  0., 2π/kx, 32
vmin, vmax, nv = -6., 6., 64

mesh_x = UniformMesh( xmin, xmax, nx )
mesh_v = UniformMesh( vmin, vmax, nv )

advection_x! = PeriodicAdvection( mesh_x, Bspline(degree))
advection_v! = PeriodicAdvection( mesh_v, Bspline(degree))

ex  = zeros(Float64, nx)
rho = zeros(Float64, nx)

fxv = zeros(Float64, (nx, nv))
fvx = zeros(Float64, (nv, nx))

x = mesh_x.points
v = mesh_v.points

fxv = (1 .+ eps .* cos.(kx .* x)) .* transpose(exp.(-.5*v.^2)) ./ sqrt(2π)

transpose!(fvx, fxv)

surface(fxv)

# +
tspan  = LinRange(0, 60, 600)
dt = 0.1

compute_charge!(rho, mesh_v, fvx)

plot(x, rho)

# +
compute_e!(ex, mesh_x, rho)

plot(x, ex)

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


