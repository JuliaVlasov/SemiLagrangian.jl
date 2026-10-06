# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Julia 1.8.5
#     language: julia
#     name: julia-1.8
# ---

using LinearAlgebra
using DoubleFloats
using SemiLagrangian
using ProgressMeter
using Plots

function printout(advd::AdvectionData, str)
    println(str)
end
printout(str) = println(str)

function trace_energy(advd::AdvectionData, t)
    if t == 0
        printout(advd, "#time\tel-energy\tkinetic-energy\tglobal-energy")
    end
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

function run_simulation(T::DataType, nbdt, sz, dt, interpall, split, tabst)
    
    epsilon = T(0.5)
    dt = T(dt)

    sp1min, sp1max, nsp1 = T(0), T(4big(pi)), sz[1]
    v1min, v1max, nv1 = -T(6.0), T(6.0), sz[3]

    mesh1_sp = UniformMesh(sp1min, sp1max, nsp1)
    mesh1_v = UniformMesh(v1min, v1max, nv1)

    sp2min, sp2max, nsp2 = T(0), T(4big(pi)), sz[2]
    v2min, v2max, nv2 = -T(6.0), T(6.0), sz[4]

    mesh2_sp = UniformMesh(sp2min, sp2max, nsp2)
    mesh2_v = UniformMesh(v2min, v2max, nv2)

    adv = Advection(
        (mesh1_sp, mesh2_sp, mesh1_v, mesh2_v),
        [interpall...],
        dt,
        tabst;
        tab_coef = split(dt),
    )

    fct_sp(x) = epsilon * cos(x / 2) + 1
    fct_v(v) = exp(-v^2 / 2) / sqrt(2T(pi))
    lgn1_sp = fct_sp.(mesh1_sp.points)
    lgn1_v = fct_v.(mesh1_v.points)
    lgn2_sp = fct_sp.(mesh2_sp.points)
    lgn2_v = fct_v.(mesh2_v.points)

    data = dotprod((lgn1_sp, lgn2_sp, lgn1_v, lgn2_v))

    pvar = getpoissonvar(adv)

    advd = AdvectionData(adv, data, pvar)

    printout(
        advd,
        "# dt=$(Float32(dt)) eps=$(Float64(epsilon)) size1_sp=$nsp1 size2_sp=$nsp2 size_v1=$nv1 size_v2=$nv2",
    )
    printout(advd, "# sp1 : from $(Float64(start(mesh1_sp))) to $(Float64(stop(mesh1_sp)))")
    printout(advd, "# sp2 : from $(Float64(start(mesh2_sp))) to $(Float64(stop(mesh2_sp)))")
    printout(advd, "# v1 : from $(Float64(start(mesh1_v))) to $(Float64(stop(mesh1_v)))")
    printout(advd, "# v2 : from $(Float64(start(mesh2_v))) to $(Float64(stop(mesh2_v)))")
    printout(advd, "# interpolation : $interpall")
    printout(advd, "# tab_coef : $split")
    printout(advd, "# tabst : $tabst")
    printout(advd, "# type=$T precision = $(precision(T))")
    printout(advd, "# monothread version")
    printout(advd, "typeof(data)=$(typeof(data)) size(data)=$(size(data))")

    maxdiff = 0
    dt = advd.adv.dt_base
    refel = getenergyall(advd)
    maxel = minel = refel
    el = Float64[]
    @showprogress 1 for i = 1:nbdt
        while advection!(advd) end
        push!(el, getenergyall(advd))
    end
    return el
end

T = Float64
nbdt = 50
sz = (32, 32, 32, 32)
dt = big"0.1"
interpall = ntuple(x -> Lagrange(19, T), 4)
split = strangsplit
tabst = [
        ([3, 4, 1, 2], 1, 1, true),
        ([4, 3, 1, 2], 1, 1, true),
        ([1, 2, 4, 3], 1, 2, true),
        ([2, 1, 3, 4], 1, 2, true),
    ]
@time result = run_simulation(T, nbdt, sz, dt, interpall, split, tabst)

t = LinRange(0, nbdt * dt, nbdt)
plot(t, result)

run_simulation(Float64, 50, (32, 32, 32, 32), big"0.1", ntuple(x -> Lagrange(19, Float64), 4), strangsplit, tabst)

run_simulation(Float64, 50, (32,32,32,32), big"0.1", ntuple(x -> Lagrange(19, Float64), 4), strangsplit, tabst)

run_simulation(BigFloat, 10000, (64,64,64,64), big"0.01", ntuple(x -> Lagrange(19, BigFloat), 4), strangsplit, tabst)

run_simulation(BigFloat, 10000, (32,32,32,32), big"0.01", ntuple(x -> Lagrange(19, BigFloat), 4), strangsplit, tabst)

T = Float64
run_simulation(T, 10000, (32,32,32,32), big"0.01", ntuple(x -> BSplineLU(27,32,T), 4), strangsplit, tabst)

@time run_simulation(T, 1000, (32,32,32,32), big"0.1", ntuple(x -> Lagrange(5, T), 4), strangsplit, tabst)

@time run_simulation(T, 30, (32,64,36,40), big"0.1", ntuple(x -> Lagrange(19, T), 4), strangsplit, tabst)

sz = (32, 32, 20, 22)

@time run_simulation(
    T,
    10,
    sz = sz,
    dt = big"0.1",
    interpall = ntuple(x -> BSplineLU(13, sz[x], T), 4),
)

@time run_simulation(T, 640, (32,32,128,128), big"0.125", ntuple(x -> Lagrange(5, T), 4), strangsplit, tabst)

c = BigFloat(2)^(1 // 3)
c1 = 1 / (2(2 - c))
c2 = (1 - c) / (2(2 - c))
d1 = 1 / (2 - c)
d2 = -c / (2 - c)
tc = [c1, d1, c2, d2, c2, d1, c1]
tc = [1, 1]
T = Float64
@time run_simulation(
    T,
    5,
    sz = (32, 32, 32, 32),
    dt = T(0.01),
    interpall = ntuple(x -> Lagrange(9, T), 4),
    split = strangsplit,
)
@time run_simulation(
    T,
    10,
    sz = (32, 32, 32, 32),
    dt = T(0.005),
    interpall = ntuple(x -> Lagrange(9, T), 4),
    split = strangsplit,
)
# +
tabst = map(
    x -> if x%2 == 1
            ([1,2,3,4], 2, x, true)
        else # x%2 == 0
            ([3,4,1,2], 2, x, true)
        end,
    1:3
)

@time run_simulation(
    T,
    30,
    sz = (32, 64, 36, 40),
    dt = big"0.1",
    interpall = ntuple(x -> Lagrange(7, T), 4),
    tabst = tabst,
)
# -
@time run_simulation(T, 10000, (64, 64, 64, 64), big"0.01", ntuple(x -> Lagrange(27, T), 4), strangsplit, tabst)
