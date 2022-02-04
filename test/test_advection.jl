
using DoubleFloats
using Random

using SemiLagrangian:
    Advection,
    AdvectionData,
    getstcoef,
    advection!,
    sizeall,
    getext,
    getdata,
    getcur_t,
    getindsplit,
    _getcurrentindice,
    getprecal,
    getitr,
    gett_split,
    nextstate!,
    UniformMesh,
    totuple,
    tovector,
    Lagrange,
    getpoissonvar,
    compute_ke,
    getinterp,
    points,
    dotprod

function initmesh(t_deb, t_end, t_size)
    t_step = (t_end - t_deb) ./ t_size
    return totuple(UniformMesh.(t_deb, t_end, t_size)), t_step
end

function test_adv(T::DataType)
    println("trace1")
    t_debsp = T.([-1, -10, -3])
    t_endsp = T.([3, 6, 5])
    t_szsp = (16, 8, 32)
    t_debv = T.([-3, -9, 1])
    t_endv = T.([1, 7, 1])
    t_szv = (4, 8, 4)
    base_dt = one(T) / 80

    Nsum = 6

    t_meshsp, t_stepsp = initmesh(t_debsp, t_endsp, t_szsp)
    t_meshv, t_stepv = initmesh(t_debv, t_endv, t_szv)

    interp = Lagrange(3, T)
    adv = Advection(
        (t_meshsp..., t_meshv...),
        map(x -> Lagrange(3, T), 1:6),
        base_dt,
        [
            ([1, 2, 3, 6, 5, 4], 1, 1, true),
            ([2, 1, 3, 4, 6, 5], 1, 1, true),
            ([3, 2, 1, 4, 5, 6], 1, 1, true), # space states
            ([4, 5, 6, 1, 2, 3], 1, 2, true),
            ([5, 4, 6, 1, 2, 3], 1, 2, true),
            ([6, 5, 4, 1, 2, 3], 1, 2, true), # velocity states
        ],
    )
    println("trace2")

    sref = (t_szsp..., t_szv...)
    @test sref == sizeall(adv)

    tab = rand(T, sizeall(adv))

    advd = AdvectionData(adv, tab, getpoissonvar(adv))

    @test compute_ke(t_meshsp, t_meshv, tab) == compute_ke(advd)

    @test advd.state_gen == 1

    @test advd.parext == getext(advd)
    @test advd.data == getdata(advd)
    @test adv.tab_coef[1] == getcur_t(adv, 1)
    @test adv.tab_coef[3] == getcur_t(adv, 1)
    @test adv.tab_coef[1] == getcur_t(adv, 2)
    @test adv.tab_coef[1] == getcur_t(adv, 3)
    @test adv.tab_coef[2] == getcur_t(adv, 4)
    @test adv.tab_coef[2] == getcur_t(adv, 5)
    @test adv.tab_coef[2] == getcur_t(adv, 6)

    @test getcur_t(advd) == advd.adv.tab_coef[1]

    println("trace3")

    t_coef = [1, 1, 1, 2, 2, 2, 3, 3, 3, 1]
    #    t_dim = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1]
    t_indice = [1, 2, 3, 4, 5, 6, 1, 2, 3, 1]
    t_v = [false, false, false, true, true, true, false, false, false, false]
    t_result = [true, true, true, true, true, true, true, true, false, true]
    # resfirst = [4, 4, 4, (4, 1, 1), (4, 1, 1), (4, 1, 1), 4, 4, 4, 4]
    # ressecond = [
    #     (:, 3, 1, 4, 1, 1),
    #     (3, :, 1, 1, 4, 1),
    #     (3, 1, :, 1, 1, 4),
    #     (4, 1, 1, :, 3, 1),
    #     (4, 1, 1, 3, :, 1),
    #     (4, 1, 1, 3, 1, :),
    #     (:, 3, 1, 4, 1, 1),
    #     (3, :, 1, 1, 4, 1),
    #     (3, 1, :, 1, 1, 4),
    #     (:, 3, 1, 4, 1, 1),
    # ]

    for i = 1:length(t_coef)
        @test getstcoef(advd) == t_coef[i]
        @test advd.state_gen == (i - 1) % 9 + 1
        @test t_indice[i] == _getcurrentindice(advd)

        @test getcur_t(advd) == adv.tab_coef[t_coef[i]]
        t = adv.t_interp
        @test t[t_indice[i]] == getinterp(advd)[1]

        # itrfirst = getitrfirst(advd)
        # (res, _) = Iterators.peel(Iterators.drop(itrfirst,3))

        # @test resfirst[i] == res

        # itrsecond = getitrsecond(advd,res)

        # (res2, _) = Iterators.peel(Iterators.drop(itrsecond,2))

        # @test ressecond[i] == res2

        #        x = t_indice[i]
        #    @time @test addcolon.(x, Iterators.product(refitr[vcat(1:(x-1),(x+1):Nsum)]...)) == getitr(advd)

        ret = nextstate!(advd)
        @test ret == t_result[i]
    end
end

# @testset "test Advection1d Float" begin

#     test_adv(Float64)

# end
# @testset "test Advection1d BigFloat" begin

#     test_adv(BigFloat)

# end
# @testset "test Advection1d Double64" begin

#     test_adv(Double64)

# end
function test_ke(T::DataType)
    t_debsp = T.([-1 // 1, -10 // 1, -3 // 1])
    t_endsp = T.([3 // 1, 6 // 1, 5 // 1])
    t_szsp = [2, 4, 8]
    t_stepsp = (t_endsp - t_debsp) ./ t_szsp
    tt_meshsp = UniformMesh.(t_debsp, t_endsp, t_szsp)
    t_meshsp = totuple(tt_meshsp)
    szsp = totuple(t_szsp)

    t_debv = T.([-3 // 1, -9 // 1, 1 // 1, -1 // 1])
    t_endv = T.([1 // 1, 7 // 1, 5 // 1, 3 // 1])
    t_szv = [4, 8, 4, 2]
    t_stepv = (t_endv - t_debv) ./ t_szv
    tt_meshv = UniformMesh.(t_debv, t_endv, t_szv)
    t_meshv = totuple(tt_meshv)
    szv = totuple(t_szv)

    fxv = if T <: Rational
        rationalize.(BigInt, rand(Float64, (szsp..., szv...)))
    else
        rand(T, (szsp..., szv...))
    end

    Nsp = length(szsp)
    Nv = length(szv)
    Nsum = Nsp + Nv
    dx = prod(step, t_meshsp)
    dv = prod(step, t_meshv)
    sum_sp = Array{T,Nv}(undef, szv)
    sum_sp .= reshape(sum(fxv; dims = ntuple(x -> x, Nsp)), szv)
    refres = (dx * dv) * sum(dotprod(points.(t_meshv)) .^ 2 .* sum_sp)

    @test refres == compute_ke(t_meshsp, t_meshv, fxv)
end

# function test_itr(T::DataType)
#     t_debsp = T.([-1, -10, -3])
#     t_endsp = T.([3, 6, 5])
#     t_szsp = (16, 8, 32)
#     t_debv = T.([-3, -9, 1])
#     t_endv = T.([1, 7, 1])
#     t_szv = (4, 8, 4)
#     base_dt = one(T) / 80

#     Nsum = 6

#     t_meshsp, t_stepsp = initmesh(t_debsp, t_endsp, t_szsp)
#     t_meshv, t_stepv = initmesh(t_debv, t_endv, t_szv)

#     interp = Lagrange(3, T)
#     adv = Advection1d(
#         t_meshsp,
#         t_meshv,
#         ntuple(x -> Lagrange(3, T), 3),
#         ntuple(x -> Lagrange(3, T), 3),
#         base_dt,
#     )

#     sref = (t_szsp..., t_szv...)

#     refitr = ntuple(x -> 1:sref[x], size(sref, 1))

#     tab = rand(T, sizeall(adv))

#     @show size(tab)

#     advd = Advection1dData(adv, tab, getpoissonvar(adv))

#     #    @show advd.t_itrfirst

#     resfirst = [
#         CartesianIndex(2, 17, 1, 2, 1),
#         CartesianIndex(2, 9, 3, 1, 1),
#         CartesianIndex(2, 1, 2, 3, 1),
#         CartesianIndex(2, 1, 5, 3, 1),
#         CartesianIndex(2, 1, 9, 5, 1),
#         CartesianIndex(2, 1, 5, 3, 1),
#         CartesianIndex(2, 17, 1, 2, 1),
#         CartesianIndex(2, 9, 3, 1, 1),
#         CartesianIndex(2, 1, 2, 3, 1),
#         CartesianIndex(2, 17, 1, 2, 1),
#     ]
#     ressecond = [
#         CartesianIndex(6, 21, 1, 3, 1),
#         CartesianIndex(6, 11, 1, 2, 1),
#         CartesianIndex(6, 5, 2, 5, 1),
#         CartesianIndex(6, 1, 6, 5, 1),
#         CartesianIndex(2, 2, 11, 1, 2),
#         CartesianIndex(6, 1, 6, 5, 1),
#         CartesianIndex(6, 21, 1, 3, 1),
#         CartesianIndex(6, 11, 1, 2, 1),
#         CartesianIndex(6, 5, 2, 5, 1),
#         CartesianIndex(6, 21, 1, 3, 1),
#     ]

#     for i = 1:length(resfirst)
#         itr = getitr(advd)

#         (res, _) = Iterators.peel(Iterators.drop(itr, 1153))

#         #        @show res

#         @test resfirst[i] == res

#         (res2, _) = Iterators.peel(Iterators.drop(itr, 2213))

#         #        @show res2

#         @test ressecond[i] == res2

#         nextstate!(advd)
#     end

# end

@testset "test compute_ke" begin
    test_ke(Rational{BigInt})
    test_ke(Float64)
    test_ke(BigFloat)
end

@testset "test advection" begin
    test_adv(Float64)
end
# @testset "test itr" begin
#     test_itr(Float64)
#     test_itr(BigFloat)
# end
