# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.5
#   kernelspec:
#     display_name: Julia 1.7.2
#     language: julia
#     name: julia-1.7
# ---

using LinearAlgebra
using DoubleFloats
using SemiLagrangian

function printout(advd::AdvectionData{T,N,timeopt}, str) where {T,N,timeopt}
    if timeopt != MPIOpt || advd.adv.mpid.ind == 1
        println(str)
    end
end
printout(str) = println(str)

function trace_energy(advd::AdvectionData{T,N,timeopt}, t) where {T,N,timeopt}
    t == 0 && printout(advd, "#time\tel-energy\tkinetic-energy\tglobal-energy")
    compute_charge!(advd)
    compute_elfield!(advd)
    elenergy = compute_ee(advd)
    kinenergy = compute_ke(advd)
    energyall = elenergy + kinenergy
    printout(
        advd,
        "$(Float32(t))\t$(Float64(elenergy))\t$(Float64(kinenergy))\t$(Float64(energyall))",
    )
    return energyall
end

function run_simulation( T::DataType, nbdt, timeop, sz, dt, interp, tab_coef)
    
    epsilon = T(0.5)
    dt = T(dt)

    spmin, spmax, nsp = T(0), T(4big(pi)), sz[1]
    vmin, vmax, nv = -T(10), T(10), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    tabst = [([1, 2], 1, 1, true), ([2, 1], 1, 2, true)]

    adv = Advection(
        (mesh_sp, mesh_v),
        [interp, interp],
        dt,
        tabst;
        tab_coef = tab_coef,
        timeopt = timeopt,
    )

    fct_sp(x) = epsilon * cos(x / 2) + 1
    fct_v(v) = exp(-v^2 / 2) / sqrt(T(2π))

    lgn_sp = fct_sp.(mesh_sp.points)
    lgn_v = fct_v.(mesh_v.points)

    data = dotprod((lgn_sp, lgn_v))

    pvar = getpoissonvar(adv)

    advd = AdvectionData(adv, data, pvar)

    printout(advd, "# dt=$(Float64(dt)) eps=$(Float64(epsilon)) size_x=$nsp size_v=$nv")
    printout(advd, "# sp : from $(Float64(start(mesh_sp))) to $(Float64(stop(mesh_sp)))")
    printout(advd, "# v : from $(Float64(start(mesh_v))) to $(Float64(stop(mesh_v)))")
    printout(advd, "# interpolation : $interp order=$(get_order(interp))")
    printout(advd, "# tab_coef : $tab_coef")
    printout(advd, "# type=$T precision = $(precision(T))")
    printout(advd, "# timeopt=$timeopt")
    if timeopt == SimpleThreadsOpt || timeopt == SplitThreadsOpt
        printout(advd, "# nb threads : $(Threads.nthreads())")
    elseif timeopt == MPIOpt
        printout(advd, "# nb process : $(adv.mpid.nb)")
    else
        printout(advd, "# monothread version")
    end
    printout(advd, "typeof(data)=$(typeof(data)) size(data)=$(size(data))")

    # advdata = Advection1dData(adv, data, pvar)

    return landau(advd, nbdt)
end

data_type = Float64
nbdt = 1000
timeopt = NoTimeOpt
sz = (64, 64)
dt = 0.1
interp = Lagrange(9, Float64)
tab_coef = [1/2, 1, 1/2]
run_simulation( data_type, nbdt, timeopt, sz, dt, interp, tab_coef)


@time run_simulation(Float64, 10000, MPIOpt, (256,256), 0.01, interp, order6split(dt))


