using DoubleFloats
using LinearAlgebra

using SemiLagrangian:
    nosplit,
    Advection,
    sizeall,
    AdvectionData,
    getdata,
    advection!,
    UniformMesh,
    AbstractInterpolation,
    Lagrange,
    B_SplineLU,
    B_SplineFFT,
    gettabmod,
    CachePrecal,
    interpolate!,
    compute_charge!,
    compute_elfield!,
    compute_ee,
    compute_ke,
    getenergy,
    dotprod,
    PoissonVar,
    getpoissonvar,
    TypePoisson,
    StdPoisson,
    StdPoisson2d,
    stdtomesh,
    TimeAlgorithm,
    NoTimeAlg,
    ABTimeAlg_ip,
    ABTimeAlg_new,
    ABTimeAlg_init,
    hamsplit_3_11,
    TimeOptimization,
    NoTimeOpt,
    SimpleThreadsOpt,
    SplitThreadsOpt



function test_poisson2d(
    sz::NTuple{2,Int},
    interp::Vector{I},
    t_max::T,
    nbdt::Int,
) where {T,I<:AbstractInterpolation{T}}
    spmin, spmax, nsp = T(0), 4T(pi), sz[1]
    vmin, vmax, nv = T(-6), T(6), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    dt = t_max / nbdt

    # dec1 = -2tan(dt/2) / step(mesh_sp) * mesh_v.points
    # dec2 = 2tan(dt/2) / step(mesh_v) * mesh_sp.points
    dec1 = zeros(T, sz)
    dec2 = zeros(T, sz)
    coef = 1 / sqrt(2T(pi))
    tabref = zeros(T, sz)
    for i = 1:sz[1], j = 1:sz[2]
        x = mesh_sp.points[i]
        y = mesh_v.points[j]
        dec1[i, j] = -dt / step(mesh_sp) * y
        tabref[i, j] = coef * exp(-0.5 * y^2) * (1 + 0.001 * cos(x / 2))
    end


    #    @show dec1, dec2

    cache = CachePrecal(interp, t_max)

    rho = zeros(T, sz[1])
    elf = zeros(T, sz[1])


    data = zeros(T, sz)

    copyto!(data, tabref)

    compute_charge!(rho, (mesh_sp,), data)
    compute_elfield!(elf, mesh_sp, rho)
    ee = compute_ee((mesh_sp,), (elf,))
    ke = compute_ke((mesh_sp,), (mesh_v,), data)

    println("0\t$ee\t$ke")


    buf = zeros(T, sz)
    diffmax = 0
    tabmod = gettabmod.(sz)
    for ind = 1:2nbdt
        dec2 = [dt / step(mesh_v) * elf[i] for i = 1:sz[1], j = 1:sz[2]]
        interpolate!(
            buf,
            data,
            ind -> (dec1[ind], dec2[ind]),
            interp,
            tabmod = tabmod,
            cache = cache,
        )
        copyto!(data, buf)

        compute_charge!(rho, (mesh_sp,), data)
        compute_elfield!(elf, mesh_sp, rho)

        ee = compute_ee((mesh_sp,), (elf,))
        ke = compute_ke((mesh_sp,), (mesh_v,), data)

        println("$(Float32(ind*dt))\t$ee\t$ke")
    end
    return diffmax
end
verif(pv::PoissonVar{T}, advd::AdvectionData) where {T} = missing
# function verif(pv::PoissonVar{T,N,Nsp,Nv,StdPoisson2dTry}, advd::AdvectionData) where{T,N,Nsp,Nv}
#     adv = advd.adv
#     sz = sizeall(adv)
#     d = zeros(T, sz)
#     coef = 1 / sqrt(2T(pi))
#     for ind in CartesianIndices(sz)
#         sp = stdtomesh(adv.t_mesh[1],pv.t_sp[ind])
#         v = stdtomesh(adv.t_mesh[2],pv.t_v[ind])
#         d[ind] = coef * exp(T(-0.5) * v^2) * (1 + T(big"0.001") * cos( sp / 2))
#     end
#     # println("t_sp=$(pv.t_sp[1:10,1:10])")
#     # println("t_v=$(pv.t_v[1:10,1:10])")
#     if !ismissing(pv.bufcur_sp)
#         # println("bufc_sp=$(pv.bufcur_sp[1:10])")
#         # println("bufc_v=$(pv.bufcur_v[1:10])")
#         println("extrema sp : $(extrema(pv.bufcur_sp))")
#         println("extrema v : $(extrema(pv.bufcur_v))")
#     end
#     res = norm(d-advd.data)
#     res2 = norm(d-advd.data, Inf)
#     println("verif : res=$res res2=$res2")
#     @show T
# end
function initdata(adv::Advection{T}) where {T}
    sz = sizeall(adv)
    coef = 1 / sqrt(2T(pi))
    tabref = zeros(T, sz)
    for i = 1:sz[1], j = 1:sz[2]
        x = adv.t_mesh[1].points[i]
        y = adv.t_mesh[2].points[j]
        tabref[i, j] = coef * exp(T(-0.5) * y^2) * (1 + T(big"0.5") * cos(x / 2))
    end
    return tabref
end


function get_init(advdref::AdvectionData{T}, nb::Int) where {T}

    dt = advdref.adv.dt_base

    tabst = [([2, 1], 1, 1, true, true), ([1, 2], 1, 2, true, false)]

    lag = Lagrange(19, T)

    adv = Advection(advdref.adv.t_mesh, [lag, lag], dt, tabst; tab_coef = hamsplit_3_11(dt))


    pvar = getpoissonvar(adv)

    advd = AdvectionData(adv, advdref.data, pvar)

    t_res = []
    for i = 1:nb
        while advection!(advd)
        end
        push!(t_res, copy(advd.data))
    end


    return t_res
end

function test_poisson2dadv(
    sz::NTuple{2,Int},
    interp::Vector{I},
    t_max::T,
    nbdt::Int,
    type::TypePoisson = StdPoisson2d,
    typeadd = 0,
    timealg::TimeAlgorithm = NoTimeAlg,
    ordalg = 0;
    timeopt::TimeOptimization=NoTimeOpt;
    flbiginit::Bool=false,
) where {T,I<:AbstractInterpolation{T}}
    spmin, spmax, nsp = T(0), 4T(pi), sz[1]
    vmin, vmax, nv = T(-9), T(9), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    dt = t_max / nbdt

    tabst = [([1, 2], 2, 1, false, false)]

    adv = Advection(
        (mesh_sp, mesh_v),
        interp,
        dt,
        tabst,
        tab_coef = nosplit(dt),
        timealg = timealg,
        ordalg = ordalg,
        timeopt=timeopt,
    )

    tabref = initdata(adv)

    data = zeros(T, sz)
    copyto!(data, tabref)

    pvar = getpoissonvar(adv, type = type, typeadd = typeadd)


    # initdatas = timealg == ABTimeAlg_init ? map(x -> tabref, 1:(ordalg-1)) : missing

    # advd = AdvectionData(adv, data, pvar; initdatas = initdatas)

    advd = AdvectionData(adv, data, pvar)

    if timealg == ABTimeAlg_init
        nb = flbiginit ? 3ordalg-1 : ordalg-1
        advd.time_cur -= dt*nb
        advd.initdatas = get_init(advd, nb)
    end

    elenergy, kinenergy, energyall = getenergy(advd)

    verif(pvar, advd)

    enmax = enmin = energyall
    @show enmax, enmin

    diff = 0
    diffprec = 0

    borne_t = t_max - dt/2

    while advd.time_cur < borne_t
        while advection!(advd)
        end
        elenergy, kinenergy, energyall = getenergy(advd)
        enmax = max(energyall, enmax)
        enmin = min(energyall, enmin)

        t = advd.time_cur
        diff = enmax - enmin
        delta = diff - diffprec
        diffprec = diff
        println("t=$(Float32(t)) diff=$diff delta=$delta")
        println(
            "$(Float32(t))\t$(Float64(elenergy))\t$(Float64(kinenergy))\t$(Float64(energyall))",
        )
        verif(pvar, advd)
    end
    @show enmax, enmin
    return (enmax - enmin), advd.data
end


function test_poisson2d2d_adv(
    sz::NTuple{4,Int},
    interp::Vector{I},
    t_max::T,
    nbdt::Int,
    split=nosplit()
) where {T,I<:AbstractInterpolation{T}}
    spmin1, spmax1, nsp1 = T(0), 4T(pi), sz[1]
    spmin2, spmax2, nsp2 = T(0), 4T(pi), sz[2]
    vmin, vmax, nv = T(-6), T(6), sz[2]
    vmin, vmax, nv = T(-6), T(6), sz[2]

    mesh_sp1 = UniformMesh(T(0), 4T(pi), sz[1])
    mesh_sp2 = UniformMesh(T(0), 4T(pi), sz[2])
    mesh_v1 = UniformMesh(T(-6), T(6), sz[3])
    mesh_v2 = UniformMesh(T(-6), T(6), sz[3])

    dt = t_max / nbdt

    tabst = [([1, 2, 3, 4], 2, 1, true, false), ([3, 4, 1, 2], 2, 2, true, true)]

    adv = Advection((mesh_sp1, mesh_sp2, mesh_v1, mesh_v2), interp, dt, tabst)
    println("trace0")

    epsilon = T(0.5)
    fct_sp(x) = epsilon * cos(x / 2) + 1
    fct_v(v) = exp(-v^2 / 2) / sqrt(2T(pi))
    lgn1_sp = fct_sp.(mesh_sp1.points)
    lgn1_v = fct_v.(mesh_v1.points)
    lgn2_sp = fct_sp.(mesh_sp2.points)
    lgn2_v = fct_v.(mesh_v2.points)
    println("trace1")
    data = dotprod((lgn1_sp, lgn2_sp, lgn1_v, lgn2_v))
    println("trace2")

    pvar = getpoissonvar(adv)

    advd = AdvectionData(adv, data, pvar)

    elenergy, kinenergy, energyall = getenergy(advd)

    enmax = enmin = energyall

    diff = 0
    diffprec = 0
    for i = 1:nbdt
        while advection!(advd)
        end
        elenergy, kinenergy, energyall = getenergy(advd)
        enmax = max(energyall, enmax)
        enmin = min(energyall, enmin)
        t = i * dt
        diff = enmax - enmin
        delta = diff - diffprec
        diffprec = diff
        println("t=$(Float32(t)) diff=$diff delta=$delta")
        println(
            "$(Float32(t))\t$(Float64(elenergy))\t$(Float64(kinenergy))\t$(Float64(energyall))",
        )
    end
    @show enmax, enmin
    return enmax - enmin
end
function test_timealg(interp::Vector{I},nbdt, timealg,ordalg) where {T,I<:AbstractInterpolation{T}}
    sz = (128,100)
   t_max = T(big"0.1")
    ret1, _ = test_poisson2dadv(
        sz,
        interp,
        t_max,
        nbdt,
        StdPoisson2d,
        0,
        timealg,ordalg)
    nbdt *= 2
    ret2, _ = test_poisson2dadv(
        sz,
        interp,
        t_max,
        nbdt,
        StdPoisson2d,
        0,
        timealg,ordalg)
    
    @show ret1,ret2, ret1/ret2

    @test 1.25*ret1/ret2 > 2 ^ ordalg
end

@testset "test poisson2d ABTimeAlg_new" begin

    T = Double64
    interp = [B_SplineLU(11, 128, T), B_SplineLU(11, 100, T)]
    @time test_timealg(interp, 5, ABTimeAlg_new, 2 )

end



@testset "test poisson2d ABTimeAlg_ip" begin
    T = Double64
    interp = [Lagrange(7,T),Lagrange(7,T)]
    @time test_timealg(interp,5,ABTimeAlg_ip,2)
    @time test_timealg(interp,5,ABTimeAlg_ip,3)
    @time test_timealg(interp,5,ABTimeAlg_ip,4)
end

@testset "test poisson2d ABTimeAlg_init" begin
    T = Double64
    interp = [Lagrange(7,T),Lagrange(7,T)]
    @time test_timealg(interp,5,ABTimeAlg_init,2)
    @time test_timealg(interp,5,ABTimeAlg_init,3)
    @time test_timealg(interp,5,ABTimeAlg_init,3)
    @time test_timealg(interp,5,ABTimeAlg_init,4)
end
@testset "test poisson2d ABTimeAlg_init 2" begin
    T = Double64
    interp = [Lagrange(9,T),Lagrange(9,T)]
    @time test_timealg(interp,5,ABTimeAlg_init,2, flbiginit=true)
    @time test_timealg(interp,5,ABTimeAlg_init,3, flbiginit=true)
    @time test_timealg(interp,5,ABTimeAlg_init,4, flbiginit=true)
    @time test_timealg(interp,10,ABTimeAlg_init,5, flbiginit=true)
end

@testset "test poisson2d thread" begin
    T = Float64
    interp = [Lagrange(5,T),Lagrange(5,T)]
    sz = (64,64)
    flag=false
    for t_alg in (NoTimeAlg, ABTimeAlg_ip, ABTimeAlg_init, ABTimeAlg_init)
        res = []
        for t_opt in (NoTimeOpt, SimpleThreadsOpt, SplitThreadsOpt)
            @show t_alg, t_opt
            ordalg = t_alg == NoTimeAlg ? 0 : 2
            _, d1 = test_poisson2dadv(
                sz,
                interp,
                T(0.01),
                5,
                StdPoisson2d,
                0,
                t_alg, ordalg,
                timeopt=t_opt,
                flbiginit=flag)
            push!(res, d1)
            if length(res) >= 2
                @test res[end-1] == res[end]
            end
            flag = t_alg == ABTimeAlg_init
        end
    end
end


