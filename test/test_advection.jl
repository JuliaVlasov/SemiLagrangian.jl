
using DoubleFloats

using SemiLagrangian:
    Advection,
    AdvectionData,
    advection!,
    sizeall,
    getext,
    getdata,
    getcur_t,
    getstate_dim,
    isvelocity,
    isvelocitystate,
    getindsplit,
    _getcurrentindice,
    getbufslgn,
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
        t_meshsp,
        t_meshv,
        ntuple(x -> Lagrange(3, T), 3),
        ntuple(x -> Lagrange(3, T), 3),
        base_dt,
    )

    sref = (t_szsp..., t_szv...)
    @test sref == sizeall(adv)


    tab = rand(T, sizeall(adv))

    advd = AdvectionData(adv, tab, getpoissonvar(adv))


    @test compute_ke(t_meshsp, t_meshv, tab) == compute_ke(advd)

    @test advd.state_coef == 1
    @test advd.state_dim == 1

    @test advd.parext == getext(advd)
    @test advd.data == getdata(advd)
    @test base_dt * adv.tab_coef[2] == getcur_t(adv, 2)
    @test 1 == getstate_dim(advd)

    @test !isvelocitystate(advd)
    @test getcur_t(advd) == base_dt * advd.adv.tab_coef[1]



    t_coef = [1, 1, 1, 2, 2, 2, 3, 3, 3, 1]
    t_dim = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1]
    t_indice = [1, 2, 3, 4, 5, 6, 1, 2, 3, 1]
    t_v = [false, false, false, true, true, true, false, false, false, false]
    t_result = [true, true, true, true, true, true, true, true, false, true]
    resfirst = [4, 4, 4, (4, 1, 1), (4, 1, 1), (4, 1, 1), 4, 4, 4, 4]
    ressecond = [
        (:, 3, 1, 4, 1, 1),
        (3, :, 1, 1, 4, 1),
        (3, 1, :, 1, 1, 4),
        (4, 1, 1, :, 3, 1),
        (4, 1, 1, 3, :, 1),
        (4, 1, 1, 3, 1, :),
        (:, 3, 1, 4, 1, 1),
        (3, :, 1, 1, 4, 1),
        (3, 1, :, 1, 1, 4),
        (:, 3, 1, 4, 1, 1),
    ]




    for i = 1:length(t_coef)
        @test advd.state_coef == t_coef[i]
        @test advd.state_dim == t_dim[i]
        @test t_indice[i] == _getcurrentindice(advd)
        @test isvelocitystate(advd) == t_v[i]
        @test advd.state_dim == getstate_dim(advd)
        @test getcur_t(advd) == base_dt * adv.tab_coef[t_coef[i]]
        @test getbufslgn(advd) == advd.t_buf[t_indice[i]]
        t = isvelocitystate(advd) ? adv.t_interp_v : adv.t_interp_sp
        @test t[t_dim[i]] == getinterp(advd)

        # itrfirst = getitrfirst(advd)
        # (res, _) = Iterators.peel(Iterators.drop(itrfirst,3))

        # @test resfirst[i] == res

        # itrsecond = getitrsecond(advd,res)

        # (res2, _) = Iterators.peel(Iterators.drop(itrsecond,2))

        # @test ressecond[i] == res2

        x = t_indice[i]
        #    @time @test addcolon.(x, Iterators.product(refitr[vcat(1:(x-1),(x+1):Nsum)]...)) == getitr(advd)

        ret = nextstate!(advd)
        @test ret == t_result[i]
    end



end

@testset "test Advection Float" begin

    test_adv(Float64)

end
@testset "test Advection BigFloat" begin

    test_adv(BigFloat)

end
@testset "test Advection Double64" begin

    test_adv(Double64)

end
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
    sum_sp .= reshape(sum(fxv, dims = ntuple(x -> x, Nsp)), szv)
    refres = (dx * dv) * sum(dotprod(points.(t_meshv)) .^ 2 .* sum_sp)

    @test refres == compute_ke(t_meshsp, t_meshv, fxv)
end

function test_itr(T::DataType)
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
        t_meshsp,
        t_meshv,
        ntuple(x -> Lagrange(3, T), 3),
        ntuple(x -> Lagrange(3, T), 3),
        base_dt,
    )

    sref = (t_szsp..., t_szv...)

    refitr = ntuple(x -> 1:sref[x], size(sref, 1))

    tab = rand(T, sizeall(adv))

    @show size(tab)

    advd = AdvectionData(adv, tab, getpoissonvar(adv))

    #    @show advd.t_itrfirst

    resfirst = [
        CartesianIndex(2, 17, 1, 2, 1),
        CartesianIndex(2, 9, 3, 1, 1),
        CartesianIndex(2, 1, 2, 3, 1),
        CartesianIndex(2, 1, 5, 3, 1),
        CartesianIndex(2, 1, 9, 5, 1),
        CartesianIndex(2, 1, 5, 3, 1),
        CartesianIndex(2, 17, 1, 2, 1),
        CartesianIndex(2, 9, 3, 1, 1),
        CartesianIndex(2, 1, 2, 3, 1),
        CartesianIndex(2, 17, 1, 2, 1),
    ]
    ressecond = [
        CartesianIndex(6, 21, 1, 3, 1),
        CartesianIndex(6, 11, 1, 2, 1),
        CartesianIndex(6, 5, 2, 5, 1),
        CartesianIndex(6, 1, 6, 5, 1),
        CartesianIndex(2, 2, 11, 1, 2),
        CartesianIndex(6, 1, 6, 5, 1),
        CartesianIndex(6, 21, 1, 3, 1),
        CartesianIndex(6, 11, 1, 2, 1),
        CartesianIndex(6, 5, 2, 5, 1),
        CartesianIndex(6, 21, 1, 3, 1),
    ]


    for i = 1:length(resfirst)
        itr = getitr(advd)

        (res, _) = Iterators.peel(Iterators.drop(itr, 1153))

        #        @show res

        @test resfirst[i] == res

        (res2, _) = Iterators.peel(Iterators.drop(itr, 2213))

        #        @show res2

        @test ressecond[i] == res2

        nextstate!(advd)
    end

end

@testset "test compute_ke" begin
    test_ke(Rational{BigInt})
    test_ke(Float64)
    test_ke(BigFloat)
end
@testset "test itr" begin
    test_itr(Float64)
    test_itr(BigFloat)
end
