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
include("../src/lagrange.jl")
include("../src/interpolation.jl")




function landau( 
    dt::T, 
    epsilon::T,
    nbdt, 
    mesh_x::UniformMesh{T}, 
    mesh_v::UniformMesh{T}, 
    interp_x::InterpolationType{T, true},
    interp_v::InterpolationType{T, true}
) where{T}
    adv_x = Advection( mesh_x, interp_x)
    adv_v = Advection( mesh_v, interp_v)

    nx = mesh_x.length
    nv = mesh_v.length

    v = mesh_v.points
    
    fxv = zeros(T, (nx, nv))
    fvx = zeros(T, (nv, nx))
    
    fct_v(v)=exp( - v^2 / 2)/sqrt(2T(pi))
    fct_x(x)=epsilon * cos(x/2) + 1
    fxv .= fct_x.(mesh_x.points) .* transpose(fct_v.(mesh_v.points))

    elf = Vector{T}(undef, nx)
    rho = Vector{T}(undef, nx)
    println("# dt=$(Float64(dt)) eps=$(Float64(epsilon)) size_x=$nx size_v=$nv")
    println("# x : from $(Float64(mesh_x.start)) to $(Float64(mesh_x.stop))")
    println("# v : from $(Float64(mesh_v.start)) to $(Float64(mesh_v.stop))")
    println("# interpolation : $(get_type(interp_x)) order=$(get_order(interp_x))")
    println("# type=$T precision = $(precision(T))")
    println("#time\tel-energy\tkinetic-energy\tglobal-energy")

    transpose!(fvx, fxv)
    compute_charge!(rho, mesh_v, fvx)
    compute_elfield!(elf, mesh_x, rho)
    elenergy = Float64(compute_ee(mesh_x, elf))
    kinenergy = Float64(compute_ke(mesh_v, mesh_x, fvx))
    energyall = elenergy + kinenergy
    println("$(Float64(0.))\t$elenergy\t$kinenergy\t$energyall")

    minall=10000000
    maxall=0
    modulo = nbdt/1000
    if modulo < 1
        modulo = 1
    end
    for i=1:nbdt
        advection!(adv_x, fxv, v, dt/2)
        transpose!(fvx, fxv)
        compute_charge!(rho, mesh_v, fvx)
        compute_elfield!(elf, mesh_x, rho)
        advection!(adv_v, fvx, elf, dt)
        transpose!(fxv, fvx)
        advection!(adv_x, fxv, v, dt/2)
        transpose!(fvx, fxv)
        compute_charge!(rho, mesh_v, fvx)
        compute_elfield!(elf, mesh_x, rho)
        elenergy = Float64(compute_ee(mesh_x, elf))
        kinenergy = Float64(compute_ke(mesh_v, mesh_x, fvx))
        energyall = elenergy + kinenergy
        if (i%modulo == 0)
            println("$(Float64(i*dt))\t$elenergy\t$kinenergy\t$energyall")
        end
        minall=min(energyall,minall)
        maxall=max(energyall,maxall)
    end
    println("diff=$(maxall-minall)")
end
    

eps    = big"0.001"
nbdt = 1000
dt = big"0.1"

xmin, xmax, nx =  big"0.", 4big(pi),  64
vmin, vmax, nv = -big"6.", big"6.", 128*8

mesh_x = UniformMesh( xmin, xmax, nx, endpoint = false, isfft=true )
mesh_v = UniformMesh( vmin, vmax, nv, endpoint = false )


interp=Lagrange(BigFloat,41)

landau(dt, eps, nbdt, mesh_x, mesh_v, interp, interp)




# ex  = zeros(Float64, nx)
# rho = zeros(Float64, nx)


# # +
# tspan  = LinRange(0, 30, 600)
# dt = 0.05

# transpose!(fvx,fxv)
# compute_charge!(rho, mesh_v, fvx)
# compute_elfield!(ex, mesh_x, rho)

# p = plot(layout=(1,2))
# scatter!(p[1,1], x, rho, label=:computed, title="rho")
# plot!(p[1,1], x, eps * cos.( kx .* x), label=:true)
# scatter!(p[1,2], x , ex, label=:computed, title="Ex")
# plot!(p[1,2], x , eps * sin.(kx .* x) / kx, label=:true)

# # +
# function simulation( fxv, tspan, mesh_x, mesh_v, interpolation )
    
#     fvx = zeros(Float64, (nv, nx))
#     transpose!(fvx, fxv)
    
#     E = Float64[]

#     @showprogress 1 for t in tspan

#         advection_x!(fxv, v, 0.5dt)
#         transpose!(fvx, fxv)
#         compute_charge!(rho, mesh_v, fvx)
#         compute_e!(ex, mesh_x, rho)
#         push!(E, 0.5 * log(sum(ex .* ex) * mesh_x.step))
#         advection_v!(fvx, ex, dt)
#         transpose!(fxv, fvx)
#         advection_x!(fxv, v, 0.5dt)

#     end
#     return E
# end

# plot(tspan, E)

# eps    = 0.001
# kx     = 0.5
# degree = 3

# xmin, xmax, nx =  0., 2π/kx, 32
# vmin, vmax, nv = -6., 6., 64

# mesh_x = UniformMesh( xmin, xmax, nx, endpoint = false )
# mesh_v = UniformMesh( vmin, vmax, nv, endpoint = false )
