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
    nosplit,
    initdata!,
    getcur_t

using DoubleFloats
using LinearAlgebra


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

    borne_t = t_max - getcur_t(advd)/2

    while advd.time_cur < borne_t
        while advection!(advd)
        end
        @show advd.time_cur
    end
    return copy(advd.data)
end
T = Double64
@time data_ref = test_quasigeostrophic((128, 128), [Lagrange(17, T), Lagrange(17, T)], T(1), 40, ABTimeAlg_ip, 6)
@time data1_10 = test_quasigeostrophic((128, 128), [Lagrange(17, T), Lagrange(17, T)], T(1), 10)
@time data1_20 = test_quasigeostrophic((128, 128), [Lagrange(17, T), Lagrange(17, T)], T(1), 20)
ret1_10 = norm(data_ref-data1_10)
ret1_20 = norm(data_ref-data1_20)
@show ret1_10,ret1_20
@show ret1_10/ret1_20

@test (1.1*ret1_10)/ret1_20 > 2

@time data2_10 = test_quasigeostrophic((128, 128), [Lagrange(17, T), Lagrange(17, T)], T(1), 10, ABTimeAlg_ip, 2)
@time data2_20 = test_quasigeostrophic((128, 128), [Lagrange(17, T), Lagrange(17, T)], T(1), 20, ABTimeAlg_ip, 2)
ret2_10 = norm(data_ref-data2_10)
ret2_20 = norm(data_ref-data2_20)
@show ret2_10,ret2_20
@show ret2_10/ret2_20
@test (1.1*ret2_10)/ret2_20 > 4


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


