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
    ABTimeAlg,
    nosplit,
    initdata!

using DoubleFloats


function test_quasigeostrophic(
    sz::NTuple{2,Int},
    interp::Vector{I},
    t_max::T,
    nbdt::Int,
    timealg::TimeAlgorithm = NoTimeAlg,
    ordalg = 0,
) where {T,I<:AbstractInterpolation{T}}
    xmin, xmax, nx = T(-5), T(5), sz[1]
    ymin, ymax, ny = T(-5), T(5), sz[2]

    mesh_x = UniformMesh(xmin, xmax, nx)
    mesh_y = UniformMesh(ymin, ymax, ny)

    dt = t_max / nbdt

    tabst = [([1, 2], 2, 1, false, false)]

    adv = Advection(
        (mesh_x, mesh_y),
        interp,
        dt,
        tabst,
        tab_coef = nosplit(dt),
        timealg = timealg,
        ordalg = ordalg,
    )

    data = zeros(T, sz)

    pvar = getgeovar(adv)


    advd = AdvectionData(adv, data, pvar)

    initdata!(pvar, advd)

    while advd.time_cur < t_max
        while advection!(advd)
        end
        @show advd.time_cur
    end
    return 0
end
T = Double64
@time res = test_quasigeostrophic((128, 128), [Lagrange(11, T), Lagrange(11, T)], T(1), 20)
@test res < 1e-10

