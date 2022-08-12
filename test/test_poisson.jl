
using DoubleFloats

using SemiLagrangian:
    UniformMesh,
    AbstractInterpolation,
    points,
    dotprod,
    compute_ee,
    compute_ke,
    PoissonConst,
    PoissonVar,
    getpoissonvar,
    initcoef!,
    compute_charge!,
    compute_elfield,
    totuple,
    Lagrange,
    Advection,
    AdvectionData,
    advection!,
    NoTimeOpt,
    SimpleThreadsOpt,
    SplitThreadsOpt,
    PrepareFftBig,
    standardsplit,
    strangsplit,
    triplejumpsplit,
    order6split,
    hamsplit_3_11,
    getenergyall

function initmesh(t_deb, t_end, t_size)
    t_step = (t_end - t_deb) ./ t_size
    return totuple(UniformMesh.(t_deb, t_end, t_size)), t_step
end

function test_poisson(T::DataType, isfft = true)
    t_debsp = T.([-1 // 1, -10 // 1, -3 // 1])
    t_endsp = T.([3 // 1, 6 // 1, 5 // 1])
    t_szsp = (8, 4, 16)
    Nsp = 3
    t_debv = T.([-3 // 1, -9 // 1, 1 // 1])
    t_endv = T.([1 // 1, 7 // 1, 5 // 1])
    t_szv = (4, 8, 4)
    base_dt = one(T) / 80
    Nv = 3

    N = Nsp + Nv

    t_meshsp, t_stepsp = initmesh(t_debsp, t_endsp, t_szsp)
    t_meshv, t_stepv = initmesh(t_debv, t_endv, t_szv)

    interp = Lagrange(3, T)
    adv = Advection(
        (t_meshsp..., t_meshv...),
        map(x -> interp, 1:6),
        base_dt,
        [([1, 2, 3, 4, 5, 6], 3, 1, false), ([4, 5, 6, 1, 2, 3], 3, 2, false)],
    )

    tab = rand(T, t_szsp..., t_szv...)

    pc = PoissonConst(adv; isfftbig = isfft)

    @test pc.v_square == dotprod(points.(adv.t_mesh[(N-Nv+1):N])) .^ 2

    # #    println("t_perms=$(pc.t_perms)")
    # @test pc.t_perms == (
    #     [1, 2, 3, 6, 5, 4],
    #     [2, 1, 3, 4, 6, 5],
    #     [3, 2, 1, 4, 5, 6], # space states
    #     [4, 5, 6, 1, 2, 3],
    #     [5, 4, 6, 1, 2, 3],
    #     [6, 5, 4, 1, 2, 3],  # velocity states
    # )

    pvar = PoissonVar(pc)

    advdata = AdvectionData(adv, tab, pvar)

    @test compute_ke(advdata) == compute_ke(totuple(t_meshsp), totuple(t_meshv), tab)

    advdata.state_gen = 2
    initcoef!(pvar, advdata)
    rhoref = zeros(T, t_szsp)
    compute_charge!(rhoref, t_meshv, tab)

    @test rhoref == pvar.rho

    pfft = if isfft
        PrepareFftBig(t_szsp, T; numdims = Nsp, dims = ntuple(x -> x, Nsp))
    else
        missing
    end

    refelfield = compute_elfield(t_meshsp, rhoref, pfft)

    #    prec = (T==BigFloat) ? 1e-70 : 1e-15

    for i = 1:Nsp
        #        @test isapprox(refelfield[i], pvar.t_elfield[i], atol=prec, rtol=prec)
        @test isapprox(refelfield[i], pvar.t_elfield[i])
    end

    #    @test isapprox(compute_ee(t_meshsp, refelfield), compute_ee(advdata), atol=prec, rtol=prec)
    @test isapprox(compute_ee(t_meshsp, refelfield), compute_ee(advdata))
end

@testset "test compute_ee" begin
    t_deb = [big"-1" // 1, -10 // 1, -3 // 1, -1 // 1]
    t_end = [big"3" // 1, 6 // 1, 5 // 1, 1 // 1]
    t_sz = [4, 4, 8, 16]
    t_step = (t_end - t_deb) ./ t_sz
    tt_mesh = UniformMesh.(t_deb, t_end, t_sz)
    t_mesh = totuple(tt_mesh)
    N = 4
    t = ntuple(x -> rationalize.(BigInt, rand(Float64, totuple(t_sz))), N)

    resref = sum(map(x -> prod(step, t_mesh) * sum(x .^ 2), t))

    @test resref == compute_ee(t_mesh, t)
end

function test_poisson_real(T::DataType, timeopt)
    dt = T(1 // 10)
    nbdt = 10
    epsilon = 0.5
    spmin, spmax, nsp = T(0), T(4big(pi)), 128
    vmin, vmax, nv = -T(10), T(10), 128

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    tabst = [([1, 2], 1, 1, true), ([2, 1], 1, 2, true)]

    interp = Lagrange(9, T)

    adv = Advection((mesh_sp, mesh_v), [interp, interp], dt, tabst; timeopt = timeopt)

    fct_sp(x) = epsilon * cos(x / 2) + 1
    fct_v(v) = exp(-v^2 / 2) / sqrt(2T(pi))

    lgn_sp = fct_sp.(mesh_sp.points)
    lgn_v = fct_v.(mesh_v.points)

    data = dotprod((lgn_sp, lgn_v))

    pvar = getpoissonvar(adv)

    advd = AdvectionData(adv, data, pvar)

    for i = 1:nbdt
        while advection!(advd)
        end
    end
    elenergy = compute_ee(advd)
    kinenergy = compute_ke(advd)

    return elenergy, kinenergy
end
function test_poisson_split(
    sz::NTuple{2,Int},
    interp::Vector{I},
    dt::T,
    nbdt::Int,
    tcoef::Vector{T},
) where {T,I<:AbstractInterpolation{T}}
    t_max = dt * nbdt
    epsilon = 0.5
    spmin, spmax, nsp = T(0), T(4big(pi)), sz[1]
    vmin, vmax, nv = -T(10), T(10), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    tabst = [([2, 1], 1, 1, true), ([1, 2], 1, 2, true)]

    adv = Advection((mesh_sp, mesh_v), interp, dt, tabst; tab_coef = tcoef)

    fct_sp(x) = epsilon * cos(x / 2) + 1
    fct_v(v) = exp(-v^2 / 2) / sqrt(2T(pi))

    lgn_sp = fct_sp.(mesh_sp.points)
    lgn_v = fct_v.(mesh_v.points)

    data = dotprod((lgn_sp, lgn_v))

    pvar = getpoissonvar(adv)

    advd = AdvectionData(adv, data, pvar)

    enall = getenergyall(advd)
    minen = maxen = enall
    for i = 1:nbdt
        while advection!(advd)
        end
        enall = getenergyall(advd)
        minen = min(enall, minen)
        maxen = max(enall, maxen)
        diff = maxen - minen
        @show advd.time_cur, diff
    end

    return maxen - minen
end
@testset "Poisson Float64" begin
    test_poisson(Float64, true)
    test_poisson(Float64, false)
end

@testset "Poisson BigFloat" begin
    test_poisson(BigFloat)
end

@testset "Poisson Double64" begin
    test_poisson(Double64)
end

@testset "Poisson Threads" begin
    res = test_poisson_real(Float64, NoTimeOpt)
    @test res == test_poisson_real(Float64, SimpleThreadsOpt)
    @test res == test_poisson_real(Float64, SplitThreadsOpt)
end

function test_split(T, nbdt, split, order)
    lag = Lagrange(19, T)
    dt = one(T) / nbdt
    sz = (64, 50)
    res1 = test_poisson_split(sz, [lag, lag], dt, nbdt, split(dt))
    nbdt *= 2
    dt = one(T) / nbdt
    res2 = test_poisson_split(sz, [lag, lag], dt, nbdt, split(dt))

    @show order, res1, res2, res1 / res2
    @test (1.05 * res1) / res2 > 2^order
end

@testset "Poisson split" begin
    @time test_split(Float64, 5, standardsplit, 1)
    @time test_split(Float64, 5, strangsplit, 2)
    @time test_split(Double64, 10, triplejumpsplit, 4)
    @time test_split(Double64, 5, order6split, 6)
    @time test_split(Double64, 5, hamsplit_3_11, 6)
end
