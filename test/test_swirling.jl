using DoubleFloats
using LinearAlgebra
#abstract type AbstractExtDataAdv end

using SemiLagrangian:
    Advection,
    sizeall,
    AdvectionData,
    getdata,
    advection!,
    UniformMesh,
    AbstractInterpolation,
    Lagrange,
    Hermite,
    BSplineLU,
    BSplineFFT,
    interpolate!,
    gettabmod,
    CachePrecal,
    AbstractExtDataAdv,
    initcoef!,
    getalpha,
    OpTuple,
    TimeAlgorithm,
    NoTimeAlg,
    ABTimeAlg_ip,
    ABTimeAlg_init

import SemiLagrangian: initcoef!, getalpha
# """

#    exact(tf, mesh1, mesh2)

#    Julia function to compute exact solution "

# ```math
# \frac{d f}{dt} +  (y \frac{df}{delta1} - x \frac{df}{delta2}) = 0 
# ```

# """
# function exact!(f, mesh1::UniformMesh{T}, mesh2::UniformMesh{T}, tf::T) where {T}
#     for (i, x) in enumerate(mesh1.points), (j, y) in enumerate(mesh2.points)
#         s, c = sincos(tf)
#         xn, yn = c * x - s * y, + s * x + c * y
#         f[i,j] = exp(-5*((xn)^2+(yn+T(6//5))^2))
#     end
# end

# function test_rotation(
#     sz::NTuple{2,Int}, 
#     interp_sp::AbstractInterpolation{T}, 
#     interp_v::AbstractInterpolation{T},
#     nbdt::Int
# ) where {T}
#     spmin, spmax, nsp =  T(-5), T(5),  sz[1]
#     vmin, vmax, nv = -T(12//2), T(9//2), sz[2]

#     mesh_sp = UniformMesh( spmin, spmax, nsp)
#     mesh_v = UniformMesh( vmin, vmax, nv)

#     dt = T(2big(pi)/nbdt)
#     adv = Advection1d((mesh_sp,), (mesh_v,), (interp_sp,), (interp_v,), dt; tab_fct=[tan, sin, tan])
# #    adv = Advection1d((mesh_sp,), (mesh_v,), (interp_sp,), (interp_v,), dt; tab_coef=[1], tab_fct=[identity])
#     sz = sizeall(adv)
#     tabref = zeros(T,sz)
#     exact!(tabref, mesh_sp, mesh_v, T(0))

#     pvar = getrotationvar(adv)

#     advdata = Advection1dData(adv, tabref, pvar)

#     diffmax = 0
#     data = getdata(advdata)
#     for ind=1:nbdt
#         while advection!(advdata) end
#         exact!(tabref, mesh_sp, mesh_v, dt*ind)
#         diff = norm(data .- tabref, Inf)
#         diffmax = max(diffmax, diff)
#         @show ind, diff
#     end
#     println("test_rotation sz=$sz interp=$interp_sp, $interp_v nbdt=$nbdt diffmax=$diffmax")
#     return diffmax

# end

function test_swirling(
    sz::NTuple{2,Int},
    ::T,
    interp::AbstractVector{I},
    nbdt::Int,
) where {T,I<:AbstractInterpolation{T}}
    spmin, spmax, nsp = T(0), T(1), sz[1]
    vmin, vmax, nv = T(0), T(1), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    t_max = T(1.5)
    dt = t_max / nbdt

    # dec1 = -2tan(dt/2) / step(mesh_sp) * mesh_v.points
    # dec2 = 2tan(dt/2) / step(mesh_v) * mesh_sp.points
    dec1 = zeros(T, sz)
    dec2 = zeros(T, sz)
    tgdt = 2tan(dt / 2)
    coef = 1 / (1 + tgdt^2 / 4)
    tabref = zeros(T, sz)
    for i = 1:sz[1], j = 1:sz[2]
        x = mesh_sp.points[i]
        y = mesh_v.points[j]
        v = [
            -cos(T(pi) * x)^2 * sin(2T(pi) * y) / step(mesh_sp),
            cos(T(pi) * y)^2 * sin(2T(pi) * x) / step(mesh_v),
        ]
        dec1[i, j], dec2[i, j] = v
        tabref[i, j] = ((1 - x)^2 + (1 - y)^2 < 0.8) ? 1 : 0
    end

    @show minimum(dec1), maximum(dec1)
    @show minimum(dec2), maximum(dec2)

    #    @show dec1, dec2

    data = zeros(T, sz)

    copyto!(data, tabref)

    diff = norm(data .- tabref)

    @show 0, diff

    tabmod = gettabmod.(sz)
    cache = CachePrecal(interp, one(T))

    buf = zeros(T, sz)
    for ind = 1:nbdt
        coef = dt * cos(T(pi) * (ind - 1) / nbdt)
        interpolate!(
            buf,
            data,
            x -> (coef * dec1[x], coef * dec2[x]),
            interp;
            tabmod = tabmod,
            cache = cache,
        )
        copyto!(data, buf)
        diff = norm(data .- tabref)
        @show ind, diff
    end
    println("test_swirling sz=$sz interp=$interp nbdt=$nbdt diff=$diff")
    return diff
end
mutable struct Swirling{T} <: AbstractExtDataAdv
    ref::Array{OpTuple{2,T},2}
    function Swirling(dec::Array{OpTuple{2,T},2}) where {T}
        return new{T}(copy(dec))
    end
end
@inbounds function initcoef!(ext::Swirling{T}, advd::AdvectionData{T}) where {T}
    coef = advd.adv.dt_base * cos(T(pi) * advd.time_cur / T(1.5))
    sz = sizeall(advd.adv)
    if ismissing(advd.bufcur)
        advd.bufcur = zeros(OpTuple{2,T}, sz)
    end
    advd.bufcur .= coef * ext.ref
    return missing
end

function test_swirling_adv(
    sz::NTuple{2,Int},
    ::T,
    interp::AbstractVector{I},
    nbdt::Int;
    timealg::TimeAlgorithm = NoTimeAlg,
    ordalg = 0,
) where {T,I<:AbstractInterpolation{T}}
    spmin, spmax, nsp = T(0), T(1), sz[1]
    vmin, vmax, nv = T(0), T(1), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    t_max = T(1.5)
    dt = t_max / nbdt

    # dec1 = -2tan(dt/2) / step(mesh_sp) * mesh_v.points
    # dec2 = 2tan(dt/2) / step(mesh_v) * mesh_sp.points
    dec = zeros(OpTuple{2,T}, sz)
    # tgdt = 2tan(dt / 2)
    # coef = 1 / (1 + tgdt^2 / 4)
    tabref = zeros(T, sz)
    for i = 1:sz[1], j = 1:sz[2]
        x = mesh_sp.points[i]
        y = mesh_v.points[j]
        v = (
            -cos(T(pi) * x)^2 * sin(2T(pi) * y) / step(mesh_sp),
            cos(T(pi) * y)^2 * sin(2T(pi) * x) / step(mesh_v),
        )
        dec[i, j] = OpTuple(v)
        tabref[i, j] = exp(-10 * ((T(1) / 4 - x)^2 + (T(1) / 2 - y)^2))
    end

    adv = Advection(
        (mesh_sp, mesh_v),
        interp,
        dt,
        [([1, 2], 2, 1, false)];
        tab_coef = [dt],
        timealg = timealg,
        ordalg = ordalg,
    )

    initdatas = missing
    if timealg == ABTimeAlg_init
        initdatas = [copy(tabref) for i = 1:(3ordalg-1)]
    end

    advd = AdvectionData(adv, tabref, Swirling(dec); initdatas = initdatas)

    if timealg == ABTimeAlg_init
        advd.time_cur -= length(initdatas) * dt
    end

    f(v, x) = v.v[x]
    @show extrema(f.(dec, (1,)))
    @show extrema(f.(dec, (2,)))

    diff = norm(getdata(advd) .- tabref)

    @show 0, diff

    for ind = 1:nbdt
        while advection!(advd)
        end
        diff = norm(getdata(advd) .- tabref)
        @show ind, diff
    end
    println("test_swirling_adv sz=$sz interp=$interp nbdt=$nbdt diff=$diff")
    return diff
end

@testset "test swirling" begin
    T = Float64
    @time @test test_swirling((100, 100), one(T), [Lagrange(9, T), Lagrange(9, T)], 50) < 15
end
@testset "test swirling_adv" begin
    T = Float64
    @time @test test_swirling_adv(
        (100, 100),
        one(T),
        [Lagrange(9, T), Lagrange(9, T)],
        50,
    ) < 5
    @time @test test_swirling_adv(
        (100, 100),
        one(T),
        [BSplineLU(9, 100, T), BSplineLU(9, 100, T)],
        50,
    ) < 5
    @time @test test_swirling_adv((100, 100), one(T), [Hermite(9, T), Hermite(9, T)], 50) <
                5

    @time @test test_swirling_adv(
        (100, 100),
        one(T),
        [Lagrange(9, T), Lagrange(9, T)],
        50;
        timealg = ABTimeAlg_init,
        ordalg = 4,
    ) < 2
    @time @test test_swirling_adv(
        (100, 100),
        one(T),
        [Hermite(9, T), Hermite(9, T)],
        50;
        timealg = ABTimeAlg_init,
        ordalg = 4,
    ) < 2
end
