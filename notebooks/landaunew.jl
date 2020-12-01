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
# using ProgressMeter
# using Plots

# +
import Base.Threads: @spawn, @sync, nthreads, threadid
include("../src/advection.jl")
include("../src/poisson.jl")
include("../src/spline.jl")
include("../src/bspline.jl")
include("../src/bsplinelu.jl")
include("../src/bsplinefft.jl")
include("../src/lagrange.jl")
include("../src/interpolation.jl")




function landau_old( 
    dt::T, 
    epsilon::T,
    nbdt, 
    tabcoef::Vector{T},
    mesh_x::UniformMesh{T}, 
    mesh_v::UniformMesh{T}, 
    interp_x::InterpolationType{T, true},
    interp_v::InterpolationType{T, true}
) where{T,ndims}
    adv_x = Advection( mesh_x, interp_x)
    adv_v = Advection( mesh_v, interp_v)

    nx = mesh_x.length
    nv = mesh_v.length

    t_coef = vcat(tabcoef, tabcoef[end-1:-1:1])
    nbcoef = size(t_coef,1)
  
    fct_v(v)=exp( - v^2 / 2)/sqrt(2T(pi))
    fct_x(x)=epsilon * cos(x/2) + 1
    lgn_x = fct_x.(mesh_x.points)
    lgn_v = fct_v.(mesh_v.points)

    v = mesh_v.points

    if ndims == 1
        nxtp = (nx)
        nvtp = (nv)
        fxv = zeros(T, (nx, nv))
        fvx = zeros(T, (nv, nx))
        fxv .= lgn_x .* reshape( lgn_v, (1,nv))
#        fxv .= fct_x.(mesh_x.points)  .* transpose(fct_v.(mesh_v.points))
        perm = [2,1]
    elseif ndims == 2
        nxtp = (nx, nx)
        nvtp = (nv, nv)
        fxv = zeros( nx, nx, nx, nv)
        fvx = zeros( nv, nv, nx, nx)
        fxv = lgn_x .* reshape(lgn_x,(1,nx)) .* reshape(lgn_v, (1,1,nv)) .* reshape(lgn_v,(1,1,1,nv))
        perm = [3, 4, 1, 2]
    else
        println("not yet implemented !!!")
    end
    elf = Array{T,ndims}(undef, nxtp)
    rho = Array{T, ndims}(undef, nxtp)
    println("# dt=$(Float64(dt)) eps=$(Float64(epsilon)) size_x=$nx size_v=$nv")
    println("# x : from $(Float64(mesh_x.start)) to $(Float64(mesh_x.stop))")
    println("# v : from $(Float64(mesh_v.start)) to $(Float64(mesh_v.stop))")
    println("# interpolation : $(get_type(interp_x)) order=$(get_order(interp_x))")
    println("# type=$T precision = $(precision(T))")
    println("#time\tel-energy\tkinetic-energy\tglobal-energy")

 #   transpose!(fvx, fxv)
    permutedims!(fvx, fxv, perm)
    compute_charge!(rho, mesh_v, fvx)
    compute_elfield!(elf, mesh_x, rho)
    elenergy = Float64(compute_ee(mesh_x, elf))
    kinenergy = Float64(compute_ke(mesh_v, mesh_x, fvx))
    energyall = elenergy + kinenergy
    println("$(Float64(0.))\t$elenergy\t$kinenergy\t$energyall")

    minall=10000000
    maxall=0
    t_coef *= dt
    fltraceend = true
    for i=1:nbdt
        # advection!(adv_x, fxv, v, dt/2)
        # transpose!(fvx, fxv)
        # compute_charge!(rho, mesh_v, fvx)
        # compute_elfield!(elf, mesh_x, rho)
        # advection!(adv_v, fvx, elf, dt)
        # transpose!(fxv, fvx)
        # advection!(adv_x, fxv, v, dt/2)
        # transpose!(fvx, fxv)
        # compute_charge!(rho, mesh_v, fvx)
        # compute_elfield!(elf, mesh_x, rho)
        # elenergy = Float64(compute_ee(mesh_x, elf))
        # kinenergy = Float64(compute_ke(mesh_v, mesh_x, fvx))
        # energyall = elenergy + kinenergy
        # if (i%modulo == 0)
        #     println("$(Float64(i*dt))\t$elenergy\t$kinenergy\t$energyall")
        # end
        # minall=min(energyall,minall)
        # maxall=max(energyall,maxall)
    
        for k=1:size(t_coef, 1)
            if k%2 == 1
                advection!(adv_x, fxv, v, t_coef[k])
                permutedims!(fvx, fxv, perm)
                if fltraceend || k != size(t_coef,1)
                    compute_charge!(rho, mesh_v, fvx)
                    compute_elfield!(elf, mesh_x, rho)
                end
            else
                advection!(adv_v, fvx, elf, t_coef[k])
                permutedims!(fxv, fvx, perm)
            end
        end
        if fltraceend
            elenergy = Float64(compute_ee(mesh_x, elf))
            kinenergy = Float64(compute_ke(mesh_v, mesh_x, fvx))
            energyall = elenergy + kinenergy
            println("$(Float64(i*dt))\t$elenergy\t$kinenergy\t$energyall")
        end
        minall=min(energyall,minall)
        maxall=max(energyall,maxall)
    end
    println("diff=$(maxall-minall)")
end

function trace_energy(advd::AdvectionData, t)

    if t == 0
        println("#time\tel-energy\tkinetic-energy\tglobal-energy")
    end
    global cl_obs
    clockbegin(cl_obs,6)
    compute_charge!(advd)
    compute_elfield!(advd)
    clockend(cl_obs,6)
    clockbegin(cl_obs,7)
    elenergy = Float64(compute_ee(advd))
    clockend(cl_obs,7)
    clockbegin(cl_obs,8)
    kinenergy = Float64(compute_ke(advd))
    clockend(cl_obs,8)
    energyall = elenergy + kinenergy
    println("$t\t$elenergy\t$kinenergy\t$energyall")
end


function landau(advd::AdvectionData, nbdt)

    global cl_obs
    clockreset(cl_obs)

    dt = advd.adv.dt_base
    trace_energy(advd, 0.0)
#    printall(cl_obs)
#    clockreset(cl_obs)
    for i=1:nbdt
        while advection!(advd)
        end
        trace_energy(advd, Float64(i*dt))
        # printall(cl_obs)
        # clockreset(cl_obs)
    end
    println("#  end")
# printall(cl_obs)
end
function landau1_1(T::DataType)
    epsilon = T(0.5)
    nbdt = 50
    dt = T(big"0.1")



    spmin, spmax, nsp =  T(0), T(4big(pi)),  64
    vmin, vmax, nv = -T(6), T(6.), 128

    mesh_sp = UniformMesh( spmin, spmax, nsp, endpoint = false)
    mesh_v = UniformMesh( vmin, vmax, nv, endpoint = false )


    interp=Lagrange(T,21)

    

    println("# dt=$(Float64(dt)) eps=$(Float64(epsilon)) size_x=$nsp size_v=$nv")
    println("# sp : from $(Float64(mesh_sp.start)) to $(Float64(mesh_sp.stop))")
    println("# v : from $(Float64(mesh_v.start)) to $(Float64(mesh_v.stop))")
    println("# interpolation : $(get_type(interp)) order=$(get_order(interp))")
    println("# type=$T precision = $(precision(T))")

    adv = Advection((mesh_sp,), (mesh_v,), (interp,), (interp,), dt)

    fct_sp(x)=epsilon * cos(x/2) + 1
    fct_v(v)=exp( - v^2 / 2)/sqrt(2T(pi))
    lgn_sp = fct_sp.(mesh_sp.points)
    lgn_v = fct_v.(mesh_v.points)

    data = dotprod((lgn_sp, lgn_v))

    pvar = getpoissonvar(adv)

    advdata = AdvectionData(adv, data, pvar; isthread=true)
    # advdata = AdvectionData(adv, data, pvar)

    landau(advdata, nbdt)
end   
function landau2_2(T::DataType, nbdt, isthread)
    epsilon = T(0.5)
    dt = T(big"0.1")

    sp1min, sp1max, nsp1 =  T(0), T(4big(pi)),  64
    v1min, v1max, nv1 = -T(6.), T(6.), 64

    mesh1_sp = UniformMesh( sp1min, sp1max, nsp1, endpoint = false)
    mesh1_v = UniformMesh( v1min, v1max, nv1, endpoint = false )

    sp2min, sp2max, nsp2 =  T(0), T(4big(pi)),  64
    v2min, v2max, nv2 = -T(6.), T(6.), 64

    mesh2_sp = UniformMesh( sp2min, sp2max, nsp2, endpoint = false)
    mesh2_v = UniformMesh( v2min, v2max, nv2, endpoint = false )

    interp=Lagrange(T,51)

    println("# dt=$(Float64(dt)) eps=$(Float64(epsilon)) size1_sp=$nsp1 size2_sp=$nsp2 size_v1=$nv1 size_v2=$nv2")
    println("# sp1 : from $(Float64(mesh1_sp.start)) to $(Float64(mesh1_sp.stop))")
    println("# sp2 : from $(Float64(mesh2_sp.start)) to $(Float64(mesh2_sp.stop))")
    println("# v1 : from $(Float64(mesh1_v.start)) to $(Float64(mesh1_v.stop))")
    println("# v2 : from $(Float64(mesh2_v.start)) to $(Float64(mesh2_v.stop))")
    println("# interpolation : $(get_type(interp)) order=$(get_order(interp))")
    println("# type=$T precision = $(precision(T))")
    if isthread
        println("# nb threads : $(Threads.nthreads())")
    else
        println("# monothread version")
    end

    adv = Advection((mesh1_sp, mesh2_sp), (mesh1_v, mesh2_v), (interp,interp,), (interp,interp,), dt)

    fct_sp(x)=epsilon * cos(x/2) + 1
    fct_v(v)=exp( - v^2 / 2)/sqrt(2T(pi))
    lgn1_sp = fct_sp.(mesh1_sp.points)
    lgn1_v = fct_v.(mesh1_v.points)
    lgn2_sp = fct_sp.(mesh2_sp.points)
    lgn2_v = fct_v.(mesh2_v.points)

    data = dotprod((lgn1_sp, lgn2_sp, lgn1_v, lgn2_v))
    
    println("typeof(data)=$(typeof(data)) size(data)=$(size(data))")

    pvar = getpoissonvar(adv)

    advdata = AdvectionData(adv, data, pvar; isthread=isthread)
    # advdata = AdvectionData(adv, data, pvar)

    landau(advdata, nbdt)
end
landau2_2(BigFloat, 1000, true)
