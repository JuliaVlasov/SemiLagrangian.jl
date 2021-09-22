using SemiLagrangian:
    getgeovar,
    AbstractInterpolation,
    Lagrange,
    UniformMesh,
    Advection,
    AdvectionData,
    advection!,
    TimeAlgorithm,
    NoTimeAlg,
    ABTimeAlg_ip,
    ABTimeAlg_init,
    nosplit,
    initdata!,
    getcur_t,
    nosplit,
    standardsplit,
    strangsplit,
    triplejumpsplit

using DoubleFloats
using LinearAlgebra
using Serialization


function test_quasigeostrophic(
    sz::NTuple{2,Int},
    interp::Vector{I},
    t_max::T,
    nbdt::Int,
    timealg::TimeAlgorithm = NoTimeAlg,
    ordalg = 0;
    split = nosplit,
    retdata::Bool = false,
    initdatas::Union{Missing,Vector{Array{T,2}}}=missing,
    ctrldatas::Union{Missing,Vector{Array{T,2}}}=missing,
) where {T,I<:AbstractInterpolation{T}}
    xmin, xmax, nx = T(0), T(1e6), sz[1]
    ymin, ymax, ny = T(0), T(1e6), sz[2]

    mesh_x = UniformMesh(xmin, xmax, nx)
    mesh_y = UniformMesh(ymin, ymax, ny)

    dt = t_max / nbdt

    if split == nosplit
        tabst = [([1, 2], 2, 1, false)]
    else
        tabst = [([1, 2], 1, 1, false), ([2, 1], 1, 2, false)]
    end

    adv = Advection(
        (mesh_x, mesh_y),
        interp,
        dt,
        tabst,
        tab_coef = split(dt),
        timealg = timealg,
        ordalg = ordalg,
    )

    data = zeros(T, sz)

    pvar = getgeovar(adv)


    advd = AdvectionData(adv, data, pvar, initdatas = initdatas)

    initdata!(pvar, advd)

    borne_t = t_max - adv.dt_base / 2

    diffmax = 0

    tabret = Array{T,2}[]

    cpt=ismissing(initdatas) ? 1 : length(initdatas)+1
    if ! ismissing(initdatas)
        @assert length(initdatas) == ordalg-1 "length(initdatas)=$(length(initdatas)) != ordalg-1  ordalg=$ordalg "
    end
    while advd.time_cur < borne_t
        while advection!(advd)
        end
        if cpt%1 == 0
            @show dt*cpt, advd.time_cur, borne_t, diffmax
        end
        if retdata
             push!(tabret, copy(advd.data))
        end
        if !ismissing(ctrldatas)
            diff = norm(advd.data-ctrldatas[cpt])
            diffmax = max(diff,diffmax)
            @show cpt, dt*cpt, advd.time_cur, borne_t, diff, diffmax
        end
        cpt += 1
        #        @show advd.time_cur, diff1, diff2
    end
    #    contourf(advd.data)

    return copy(advd.data)
    # println("trace1")
    # if retdata
    #     @show typeof(tabret)
    #     ret = ismissing(ctrldatas) ? tabret : (tabret, diffmax)
    #     @show typeof(ret)
    #     return ret
    # else
    #     return diffmax
    # end
end


function test_order(T, nbdt, split, str)
    m1ref, m2ref, dref = test_quasigeostrophic(
        (128, 128),
        [Lagrange(9, T), Lagrange(9, T)],
        T(10000),
        nbdt * 4,
        NoTimeAlg,
        0,
        split,
    )
    m11, m21, d1 = test_quasigeostrophic(
        (128, 128),
        [Lagrange(9, T), Lagrange(9, T)],
        T(10000),
        nbdt,
        NoTimeAlg,
        0,
        split,
    )
    m12, m22, d2 = test_quasigeostrophic(
        (128, 128),
        [Lagrange(9, T), Lagrange(9, T)],
        T(10000),
        nbdt * 2,
        NoTimeAlg,
        0,
        split,
    )
    ret1 = norm(dref - d1)
    ret2 = norm(dref - d2)

    retm1 = m11 / m12
    retm2 = m21 / m22

    @show str, ret1, ret2, ret1 / ret2
    @show retm1, retm2
end
function test_orderno(T, sz, interps, nbdt, timealg, ordalg, str)


    t_max = T(10000)


    data4 = test_quasigeostrophic(
        sz,
        interps,
        t_max,
        nbdt*4,
        timealg,
        ordalg,
    )
    data1 = test_quasigeostrophic(
        sz,
        interps,
        t_max,
        nbdt,
        timealg,
        ordalg,
    )

    data2 = test_quasigeostrophic(
        sz,
        interps,
        t_max,
        nbdt * 2,
        timealg,
        ordalg,
    )
    ret1 = norm(data4-data1)
    ret2 = norm(data4-data2)

    @show str, ret1, ret2, ret1 / ret2

    ord = ordalg == 0 ? 1 : ordalg

    @test ret1*1.2/ret2 > 2^ord
end
T = Double64
# _,_, data = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(100000), 100)
## old

# @time _,_,data_ref = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 320, ABTimeAlg_ip, 2)
sz = (128,128)
interps=[Lagrange(9,T),Lagrange(9,T)]
nbdt = 20

nbprec = 32
@show nbprec, nbdt


# @time dref = test_quasigeostrophic(
#     sz, interps,
#     T(100000),
#     2*nbdt * nbprec,
#     ABTimeAlg_ip,
#     2,
#     retdata = true,
# )

# dref = deserialize("ABTimeAlg_ip2_32.sqg")

# dref128 = deserialize("NoTimeAlg_128.sqg")

# diffip2_No = norm(norm.(dref .- dref128[4:4:end]), Inf)

# @show diffip2_No


# dr = dref[nbprec:nbprec:end]
# dr2 = dref128[128:128:end]

# # @time test_orderno(T, sz, interps, 16*nbdt, NoTimeAlg, 0, "NoTimeAlg dref", ctrldatas=dref[2:2:end])

# tabres = Array{T,2}[]
# for i=0:6
#     coef = 2 ^ i
#     @time ref = test_quasigeostrophic(
#         sz,
#         interps,
#         T(100000),
#         nbdt*coef,
#         ABTimeAlg_ip,
#         2,
#         retdata=true,
#     )
#     @show coef, typeof(ref)
#     push!(tabres, ref[end])
# end
# tabresip = deserialize("res6ip.sgq")

# tabres = deserialize("res6.sgq")
# for i = 1:7
#     for j=1:7
#         print(" $(norm(tabres[i]-tabresip[j]))")
#     end
#     println("")
# end
# for i = 1:7
#     for j=1:7
#         print(" $(log2(norm(tabres[i]-tabresip[j])))")
#     end
#     println("")
# end
# for i = 1:7
#     for j=1:7
#         print(" $(norm(tabresip[i]-tabresip[j]))")
#     end
#     println("")
# end
# for i = 1:7
#     for j=1:7
#         print(" $(log2(norm(tabresip[i]-tabresip[j])))")
#     end
#     println("")
# end
nbdt = 10

@time test_orderno(T, sz, interps, nbdt, NoTimeAlg, 0, "NoTimeAlg")
@time test_orderno(T, sz, interps,nbdt, ABTimeAlg_ip, 2, "ABTimeAlg_ip")
@time test_orderno(T, sz, interps,nbdt, ABTimeAlg_ip, 3, "ABTimeAlg_ip")

# @time test_order(T, 10, standardsplit, "standardsplit")

# @time test_order(T, 10, strangsplit, "strangsplit")
# @time test_order(T, 10, triplejumpsplit, "triplejumpsplit")



# @time test_init(T, 20, ABTimeAlg_init, 2)

# @time _,_,data1_10 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 10)
# @time _,_,data1_20 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 20)
# ret1_10 = norm(data_ref-data1_10)
# ret1_20 = norm(data_ref-data1_20)
# @show ret1_10,ret1_20
# @show ret1_10/ret1_20

# @test (1.1*ret1_10)/ret1_20 > 2

# @time _,_,data2_10 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 10, ABTimeAlg_ip, 2)
# @time _,_,data2_20 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 20, ABTimeAlg_ip, 2)
# ret2_10 = norm(data_ref-data2_10)
# ret2_20 = norm(data_ref-data2_20)
# @show ret2_10,ret2_20
# @show ret2_10/ret2_20
# @test (1.1*ret2_10)/ret2_20 > 4

# @time _,_,data3_10 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 10, ABTimeAlg_ip, 3)
# @time _,_,data3_20 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 20, ABTimeAlg_ip, 3)
# ret3_10 = norm(data_ref-data3_10)
# ret3_20 = norm(data_ref-data3_20)
# @show ret3_10,ret3_20
# @show ret3_10/ret3_20
# @test (1.1*ret3_10)/ret3_20 > 8
# ## end old

# T = Double64
# @time d1_10,d2_10 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 10)
# @time d1_20,d2_20 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 20)

# @show d1_10,d1_20, d1_10/d1_20
# @show d2_10, d2_20, d2_10/d2_20

# @test (1.1*d1_10)/d1_20 > 2
# @test (1.1*d2_10)/d2_20 > 2

# @time d1_2_10, d2_2_10 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 10, ABTimeAlg_ip, 2)
# @time d1_2_20, d2_2_20 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 20, ABTimeAlg_ip, 2)
# @show d1_2_10,d1_2_20, d1_2_10/d1_2_20
# @show d2_2_10, d2_2_20, d2_2_10/d2_2_20

# @test (1.1*d1_2_10)/d1_2_20 > 4
# @test (1.1*d2_2_10)/d2_2_20 > 4

# @time d1_3_10, d2_3_10 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 10, ABTimeAlg_ip, 3)
# @time d1_3_20, d2_3_20 = test_quasigeostrophic((128, 128), [Lagrange(9, T), Lagrange(9, T)], T(1), 20, ABTimeAlg_ip, 3)
# @show d1_3_10,d1_3_20, d1_3_10/d1_3_20
# @show d2_3_10, d2_3_20, d2_3_10/d2_3_20

# @test (1.1*d1_3_10)/d1_3_20 > 8
# @test (1.1*d2_3_10)/d2_3_20 > 8

# @time data3_10 = test_quasigeostrophic((128, 128), [Lagrange(17, T), Lagrange(17, T)], T(1), 10, ABTimeAlg_ip, 3)
# @time data3_20 = test_quasigeostrophic((128, 128), [Lagrange(17, T), Lagrange(17, T)], T(1), 20, ABTimeAlg_ip, 3)
# ret3_10 = norm(data_ref-data3_10)
# ret3_20 = norm(data_ref-data3_20)
# @show ret3_10,ret3_20
# @show ret3_10/ret3_20
# @test (1.1*ret3_10)/ret3_20 > 8

# @time data4_10 = test_quasigeostrophic((128, 128), [Lagrange(17, T), Lagrange(17, T)], T(1), 10, ABTimeAlg_ip, 4)
# @time data4_20 = test_quasigeostrophic((128, 128), [Lagrange(17, T), Lagrange(17, T)], T(1), 20, ABTimeAlg_ip, 4)
# ret4_10 = norm(data_ref-data4_10)
# ret4_20 = norm(data_ref-data4_20)
# @show ret4_10,ret4_20
# @show ret4_10/ret4_20
# @test (1.1*ret4_10)/ret4_20 > 16

# @time data5_10 = test_quasigeostrophic((128, 128), [Lagrange(17, T), Lagrange(17, T)], T(1), 10, ABTimeAlg_ip, 5)
# @time data5_20 = test_quasigeostrophic((128, 128), [Lagrange(17, T), Lagrange(17, T)], T(1), 20, ABTimeAlg_ip, 5)
# ret5_10 = norm(data_ref-data5_10)
# ret5_20 = norm(data_ref-data5_20)
# @show ret5_10,ret5_20
# @show ret5_10/ret5_20




