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
using DoubleFloats
# using SemiLagMPI

using SemiLagrangian

# function landau_old(
#     dt::T,
#     epsilon::T,
#     nbdt,
#     tabcoef::Vector{T},
#     mesh_x::UniformMesh{T},
#     mesh_v::UniformMesh{T},
#     interp_x::AbstractInterpolation{T},
#     interp_v::AbstractInterpolation{T},
# ) where {T,ndims}
#     adv_x = Advection1d(mesh_x, interp_x)
#     adv_v = Advection1d(mesh_v, interp_v)

#     nx = mesh_x.length
#     nv = mesh_v.length

#     t_coef = vcat(tabcoef, tabcoef[end-1:-1:1])
#     nbcoef = size(t_coef, 1)

#     fct_v(v) = exp(-v^2 / 2) / sqrt(2T(pi))
#     fct_x(x) = epsilon * cos(x / 2) + 1
#     lgn_x = fct_x.(mesh_x.points)
#     lgn_v = fct_v.(mesh_v.points)

#     v = mesh_v.points

#     if ndims == 1
#         nxtp = (nx)
#         nvtp = (nv)
#         fxv = zeros(T, (nx, nv))
#         fvx = zeros(T, (nv, nx))
#         fxv .= lgn_x .* reshape(lgn_v, (1, nv))
#         #        fxv .= fct_x.(mesh_x.points)  .* transpose(fct_v.(mesh_v.points))
#         perm = [2, 1]
#     elseif ndims == 2
#         nxtp = (nx, nx)
#         nvtp = (nv, nv)
#         fxv = zeros(nx, nx, nx, nv)
#         fvx = zeros(nv, nv, nx, nx)
#         fxv =
#             lgn_x .* reshape(lgn_x, (1, nx)) .* reshape(lgn_v, (1, 1, nv)) .*
#             reshape(lgn_v, (1, 1, 1, nv))
#         perm = [3, 4, 1, 2]
#     else
#         println("not yet implemented !!!")
#     end
#     elf = Array{T,ndims}(undef, nxtp)
#     rho = Array{T,ndims}(undef, nxtp)
#     printout(advd, "# dt=$(Float64(dt)) eps=$(Float64(epsilon)) size_x=$nx size_v=$nv")
#     printout(advd, "# x : from $(Float64(mesh_x.start)) to $(Float64(mesh_x.stop))")
#     printout(advd, "# v : from $(Float64(mesh_v.start)) to $(Float64(mesh_v.stop))")
#     printout(advd, "# interpolation : $(get_type(interp_x)) order=$(get_order(interp_x))")
#     printout(advd, "# type=$T precision = $(precision(T))")
#     printout(advd, "#time\tel-energy\tkinetic-energy\tglobal-energy")

#     #   transpose!(fvx, fxv)
#     permutedims!(fvx, fxv, perm)
#     compute_charge!(rho, mesh_v, fvx)
#     compute_elfield!(elf, mesh_x, rho)
#     elenergy = Float64(compute_ee(mesh_x, elf))
#     kinenergy = Float64(compute_ke(mesh_v, mesh_x, fvx))
#     energyall = elenergy + kinenergy
#     printout(advd, "$(Float64(0.))\t$elenergy\t$kinenergy\t$energyall")

#     minall = 10000000
#     maxall = 0
#     t_coef *= dt
#     fltraceend = true
#     for i = 1:nbdt
#         # advection!(adv_x, fxv, v, dt/2)
#         # transpose!(fvx, fxv)
#         # compute_charge!(rho, mesh_v, fvx)
#         # compute_elfield!(elf, mesh_x, rho)
#         # advection!(adv_v, fvx, elf, dt)
#         # transpose!(fxv, fvx)
#         # advection!(adv_x, fxv, v, dt/2)
#         # transpose!(fvx, fxv)
#         # compute_charge!(rho, mesh_v, fvx)
#         # compute_elfield!(elf, mesh_x, rho)
#         # elenergy = Float64(compute_ee(mesh_x, elf))
#         # kinenergy = Float64(compute_ke(mesh_v, mesh_x, fvx))
#         # energyall = elenergy + kinenergy
#         # if (i%modulo == 0)
#         #     println("$(Float64(i*dt))\t$elenergy\t$kinenergy\t$energyall")
#         # end
#         # minall=min(energyall,minall)
#         # maxall=max(energyall,maxall)

#         for k = 1:size(t_coef, 1)
#             if k % 2 == 1
#                 advection!(adv_x, fxv, v, t_coef[k])
#                 permutedims!(fvx, fxv, perm)
#                 if fltraceend || k != size(t_coef, 1)
#                     compute_charge!(rho, mesh_v, fvx)
#                     compute_elfield!(elf, mesh_x, rho)
#                 end
#             else
#                 advection!(adv_v, fvx, elf, t_coef[k])
#                 permutedims!(fxv, fvx, perm)
#             end
#         end
#         if fltraceend
#             elenergy = Float64(compute_ee(mesh_x, elf))
#             kinenergy = Float64(compute_ke(mesh_v, mesh_x, fvx))
#             energyall = elenergy + kinenergy
#             printout(advd, "$(Float64(i*dt))\t$elenergy\t$kinenergy\t$energyall")
#         end
#         minall = min(energyall, minall)
#         maxall = max(energyall, maxall)
#     end
#     printout(advd, "diff=$(maxall-minall)")
# end
function printout(
    advd::AdvectionData{T,N,timeopt},
    str,
) where {T,N,timeopt}
    if timeopt != MPIOpt || advd.adv.mpid.ind == 1
        println(str)
    end
end
printout(str) = println(str)

function trace_energy(
    advd::AdvectionData{T,N,timeopt},
    t,
) where {T,N,timeopt}

    if t == 0
        printout(advd, "#time\tel-energy\tkinetic-energy\tglobal-energy")
    end
    # global cl_obs
    # clockbegin(cl_obs,6)
    compute_charge!(advd)
    compute_elfield!(advd)
    # clockend(cl_obs,6)
    # clockbegin(cl_obs,7)
    elenergy = compute_ee(advd)
    # clockend(cl_obs,7)
    # clockbegin(cl_obs,8)
    kinenergy = compute_ke(advd)
    # clockend(cl_obs,8)
    energyall = elenergy + kinenergy
    printout(
        advd,
        "$(Float32(t))\t$(Float64(elenergy))\t$(Float64(kinenergy))\t$(Float64(energyall))",
    )
end


function landau(advd::AdvectionData, nbdt)

    # global cl_obs
    # clockreset(cl_obs)

    dt = advd.adv.dt_base
    trace_energy(advd, 0.0)
    #    printall(cl_obs)
    #    clockreset(cl_obs)
    for i = 1:nbdt
        while advection!(advd)
        end
        trace_energy(advd, Float64(i * dt))
        # printall(cl_obs)
        # clockreset(cl_obs)
    end
    println("#  end")
    # printall(cl_obs)
end
function landau1_1(
    T::DataType,
    nbdt,
    timeopt;
    sz = (64, 64),
    dt = big"0.1",
    interp = Lagrange(21, T),
    tab_coef=[1//2, 1, 1//2]
)
    epsilon = T(0.001)
    dt = T(dt)

    spmin, spmax, nsp = T(0), T(4big(pi)), sz[1]
    vmin, vmax, nv = -T(5), T(5), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    tabst = map(x -> ((x%2==1) ? [1,2] : [2,1], 1, x, true),1:size(tab_coef,1))

    adv = Advection(
        (mesh_sp,mesh_v,), 
        [interp,interp], 
        dt,
        tabst,
        tab_coef=tab_coef, 
        timeopt = timeopt)

    fct_sp(x) = epsilon * cos(x / 2) + 1
    fct_v(v) = exp(-v^2 / 2) / sqrt(2T(pi))

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

    landau(advd, nbdt)
end
function landau2_2(
    T::DataType,
    nbdt,
    timeopt;
    sz = (32, 32, 32, 32),
    dt = big"0.1",
    interpall = ntuple(x -> Lagrange(19, T), 4),
    tab_coef=[1//2, 1, 1//2],
    tabst = map( 
        x -> if x%4 == 1
                ([1,2,4,3], 1, div(x+1, 2), true)
            elseif x%4 == 2
                ([2,1,3,4], 1,div(x+1, 2), true)
            elseif x%4 == 3
                ([3,4,1,2], 1,div(x+1,2), true)
            else # x%4 == 0
                ([4,3,1,2], 1,div(x+1,2), true)
            end, 
        1:6
)
)
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
        tabst,
        tab_coef=tab_coef, 
        timeopt = timeopt,
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
    # advdata = Advection1dData(adv, data, pvar)
    printout(
        advd,
        "# dt=$(Float32(dt)) eps=$(Float64(epsilon)) size1_sp=$nsp1 size2_sp=$nsp2 size_v1=$nv1 size_v2=$nv2",
    )
    printout(
        advd,
        "# sp1 : from $(Float64(start(mesh1_sp))) to $(Float64(stop(mesh1_sp)))",
    )
    printout(
        advd,
        "# sp2 : from $(Float64(start(mesh2_sp))) to $(Float64(stop(mesh2_sp)))",
    )
    printout(
        advd,
        "# v1 : from $(Float64(start(mesh1_v))) to $(Float64(stop(mesh1_v)))",
    )
    printout(
        advd,
        "# v2 : from $(Float64(start(mesh2_v))) to $(Float64(stop(mesh2_v)))",
    )
    printout(advd, "# interpolation : $interpall")
    printout(advd, "# tab_coef : $tab_coef")
    printout(advd, "# tabst : $tabst")
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

    landau(advd, nbdt)
end
# landau2_2(Float64, 50, NoTimeOpt)
# landau2_2(Float64, 50, SimpleThreadsOpt)
# landau2_2(Float64, 50, SplitThreadsOpt)
# landau2_2(Float64, 50, MPIOpt, sz=(32,32,32,32))
# landau2_2(BigFloat, 10000, MPIOpt, sz=(64,64,64,64), dt=big"0.01")
# landau2_2(BigFloat, 10000, MPIOpt, sz=(32,32,32,32), dt=big"0.01")
T = Float64
# landau2_2(T, 10000, NoTimeOpt, sz=(32,32,32,32), dt=big"0.01", interp=B_SplineLU(27,32,T))
# @time landau2_2(T, 1000, NoTimeOpt, sz=(32,32,32,32), dt=big"0.1", interp=Lagrange(5, T))
# @time landau2_2(T, 30, NoTimeOpt, sz=(32,64,36,40), dt=big"0.1")
# sz = (32, 32, 20, 22)
# @time landau2_2(
#     T,
#     10,
#     MPIOpt,
#     sz = sz,
#     dt = big"0.1",
#     interpall = ntuple(x -> B_SplineLU(13, sz[x], T), 4),
# )
## @time landau2_2(T, 640, NoTimeOpt, sz=(32,32,128,128), dt=big"0.125", interp=Lagrange(5, T))

c=BigFloat(2)^(1//3)
c1 = 1/(2(2-c))
c2 = (1-c)/(2(2-c))
d1 = 1/(2-c)
d2 = -c/(2-c)
tc = [c1, d1, c2, d2, c2, d1, c1]
tc = [1, 1]
@time landau2_2(T, 30, NoTimeOpt, sz=(32,64,36,40), dt=big"0.1", interpall=ntuple(x->LagrangeInt(7,T),4),
 tabst = map( 
    x -> if x%2 != 1
            ([1,2,3,4], 2, 1, true)
        else # x%2 == 0
            ([3,4,1,2], 2, 2, true)
        end, 
    1:2
),
tab_coef=tc   
)
# @time landau1_1(T, 1000, NoTimeOpt, sz=(512,512), dt=big"0.01", tab_coef=tc)
# @time landau2_2(T, 30, NoTimeOpt, sz=(32,64,36,40), dt=big"0.1", interpall=ntuple(x->LagrangeInt(7,T),4))
# @time landau2_2(T, 30, NoTimeOpt, sz=(32,64,36,40), dt=big"0.1", interpall=ntuple(x->LagrangeInt(7,T),4),
#  tabst = map( 
#     x -> if x%2 == 1
#             ([1,2,4,3], 2, x, true)
#         else # x%2 == 0
#             ([3,4,1,2], 2, x, true)
#         end, 
#     1:3
# )   
# )
# landau1_1(T, 50, NoTimeOpt, sz=(64,128))
# landau2_2(T, 10000, MPIOpt, sz=(64,64,64,64), dt=big"0.01", interp=Lagrange(27, T))
