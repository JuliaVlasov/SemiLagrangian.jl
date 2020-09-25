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
include("../src/advection.jl")
include("../src/spline.jl")
include("../src/bspline.jl")
include("../src/bsplinelu.jl")
include("../src/bsplinefft.jl")




function landau( 
    dt::T, 
    epsilon::T, 
    mesh_x::UniformMesh{T}, 
    mesh_v::UniformMesh{T}, 
    interp_x::InterpolationType{T, true},
    interp_y::InterpolationType{T, true}
) where{T}
    advection_x! = Advection( mesh_x, interp_x)
    advection_v! = Advection( mesh_v, interp_y)

    nx = mesh_x.length
    ny = mesh_y.length
    
    fxv = zeros(T, (nx, nv))
    fvx = zeros(T, (nv, nx))
    
    fct_v(v)=exp( - v^2 / 2)/sqrt(2T(pi))
    fct_x(x)=epsilon * cos(x) + 1
    fxv .= fct_x.(mesh_x.points) .* transpose(fct_v.(mesh_v.points))

    ex = Vector{T}(nx,undef)
    rho = Vector{T}(nx,undef)


    for i=1:100
        advection_x!(fxv, v, 0.5dt)
        transpose!(fvx, fxv)
        compute_charge!(rho, mesh_v, fvx)
        compute_e!(ex, mesh_x, rho)
        advection_v!(fvx, ex, dt)
        transpose!(fxv, fvx)
        advection_x!(fxv, v, 0.5dt)
    end
end
    

    eps    = 0.001
kx     = 0.5
degree = 3

xmin, xmax, nx =  0., 2π/kx, 32
vmin, vmax, nv = -6., 6., 64

mesh_x = UniformMesh( xmin, xmax, nx, endpoint = false )
mesh_v = UniformMesh( vmin, vmax, nv, endpoint = false )


ex  = zeros(Float64, nx)
rho = zeros(Float64, nx)


# +
tspan  = LinRange(0, 30, 600)
dt = 0.05

transpose!(fvx,fxv)
compute_charge!(rho, mesh_v, fvx)
compute_elfield!(ex, mesh_x, rho)

p = plot(layout=(1,2))
scatter!(p[1,1], x, rho, label=:computed, title="rho")
plot!(p[1,1], x, eps * cos.( kx .* x), label=:true)
scatter!(p[1,2], x , ex, label=:computed, title="Ex")
plot!(p[1,2], x , eps * sin.(kx .* x) / kx, label=:true)

# +
function simulation( fxv, tspan, mesh_x, mesh_v, interpolation )
    
    fvx = zeros(Float64, (nv, nx))
    transpose!(fvx, fxv)
    
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
    return E
end

plot(tspan, E)

eps    = 0.001
kx     = 0.5
degree = 3

xmin, xmax, nx =  0., 2π/kx, 32
vmin, vmax, nv = -6., 6., 64

mesh_x = UniformMesh( xmin, xmax, nx, endpoint = false )
mesh_v = UniformMesh( vmin, vmax, nv, endpoint = false )
