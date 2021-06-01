using DoubleFloats
using LinearAlgebra

using SemiLagrangian:
    Advection,
    standardsplit,
    strangsplit,
    triplejumpsplit,
    sizeall,
    AdvectionData,
    getdata,
    advection!,
    UniformMesh,
    gettranslationvar,
    AbstractInterpolation,
    interpolate!,
    Lagrange,
    B_SplineLU,
    B_SplineFFT,
    Hermite
"""

   exact(tf, mesh1, mesh2)

   Julia function to compute exact solution "

```math
\frac{d f}{dt} +  (y \frac{df}{delta1} - x \frac{df}{delta2}) = 0 
```

"""
function exact!(f, decall::NTuple{Nsum,T}, t::T) where {T,Nsum}
    for ind in CartesianIndices(f)
        f[ind] = exp(
            -sum([
                sin(2T(pi) * (ind.I[i] + decall[i] * t) / size(f, i) + T(big"0.56") * i^2) for i = 1:Nsum
            ]),
        )
    end
end

function test_translation(
    sz::NTuple{2,Int},
    interp_sp::AbstractInterpolation{T},
    interp_v::AbstractInterpolation{T},
    dt::T;
    tc::Vector{T}=strangsplit(dt),
) where {T}
@time begin
    spmin, spmax, nsp = T(-5), T(5), sz[1]
    vmin, vmax, nv = -T(5.55), T(5), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)
end
println("trace1")
@time begin
    nbdt = Int(round(T(1)/dt))
    adv = Advection(
        (mesh_sp, mesh_v),
        [interp_sp, interp_v],
        dt,
        [([1, 2], 1, 1, true, false), ([2, 1], 1, 2, true, false)],
        tab_coef= tc
    )

    tabref = zeros(T, sz)

    v1 = T(big"0.83545655467782872872782870029282982828737872878776717190927267611111")
    v2 = T(big"0.945678101929276765616767176761671771766717828781828998101092981877817176")
end
println("trace2")
@time begin

    exact!(tabref, (v1, v2), T(0))

    pvar = gettranslationvar((v1, v2))

    advdata = AdvectionData(adv, tabref, pvar)

    diffmax = 0
    data = getdata(advdata)
end
println("trace3")
    for ind = 1:nbdt
        data = getdata(advdata)
        while advection!(advdata)
        end
        #        GC.gc()
        exact!(tabref, (v1, v2), T(ind)*dt)
        data = getdata(advdata)
        diff = norm(data .- tabref, Inf)
        @show v1, v2, ind, diff
        diffmax = max(diffmax, diff)
    end
    println(
        "length(tab_coef)=$(length(tc)) test_translation sz=$sz interp=$interp_sp, $interp_v nbdt=$nbdt diffmax=$diffmax",
    )
    return diffmax

end
function test_translation(
    sz::NTuple{N,Int},
    interp_t::Vector{I},
    nbdt::Int,
) where {T,N,I<:AbstractInterpolation{T}}

    dt = T(1)

    vall = ntuple(x -> T(10rand(BigFloat) - 5), N)

    vall = T(1) / 10 .* vall

    @show vall

    tabref = zeros(T, sz)
    exact!(tabref, vall, T(0))

    data = copy(tabref)

    buf = zeros(T, sz)

    # dec1 = fill(v1, sz)
    # dec2 = fill(v2, sz)


    diffmax = 0
    for ind = 1:nbdt
        interpolate!(buf, data, (i, _) -> vall[i], interp_t)
        #        interpolate!(buf, data, (dec1, dec2), interp)
        copyto!(data, buf)
        exact!(tabref, vall, T(ind))
        diff = norm(data .- tabref, Inf)
        diffmax = max(diffmax, diff)
        @show ind, diff
    end
    println("test_translation sz=$sz interp=$interp_t nbdt=$nbdt diffmax=$diffmax")
    return diffmax

end


@testset "test translation" begin
    T = Float64
    dt = T(1)/11
    # @time @test test_translation((150, 220), [Lagrange(9, T), Lagrange(9, T)], 11) < 1e-6
    # @time @test test_translation((50, 30, 40), map(x -> Lagrange(5 + 2x, T), 1:3), 11) <
    #             1e-4
    # @time @test test_translation((100, 120), [Lagrange(5, T), Lagrange(5, T)], 11) < 1e-5
    # @time @test test_translation(
    #     (128, 64),
    #     [B_SplineLU(9, 128, T), B_SplineLU(9, 64, T)],
    #     11,
    # ) < 1e-6
    @time @test test_translation(
        (200, 300),
        Lagrange(5, T),
        Lagrange(5, T),
        dt,
    ) < 1e-3
    @time @test test_translation(
        (200, 300),
        Lagrange(5, T),
        Lagrange(5, T),
        dt, tc=triplejumpsplit(dt)
    ) < 1e-3
    @time @test test_translation(
        (200, 300),
        Lagrange(5, T),
        Lagrange(5, T),
        dt, tc=standardsplit(dt)
    ) < 1e-3
    @time @test test_translation(
        (200, 300),
        Hermite(5, T),
        Hermite(5, T),
        dt,
    ) < 1e-3
    @time @test test_translation(
        (200, 300),
        Hermite(5, T),
        Hermite(5, T),
        dt, triplejumpsplit(dt)
    ) < 1e-3
    @time @test test_translation(
        (200, 300),
        Hermite(5, T),
        Hermite(5, T),
        dt,standardsplit(dt)
    ) < 1e-3
    @time @test test_translation(
        (128, 64),
        B_SplineLU(5, 128, T),
        B_SplineLU(5, 64, T),
        dt,
    ) < 1e-6
    @time @test test_translation(
        (128, 64),
        B_SplineLU(5, 128, T),
        B_SplineLU(5, 64, T),
        dt, triplejumpsplit(dt)
    ) < 1e-6
    @time @test test_translation(
        (128, 64),
        B_SplineFFT(5, 128, T),
        B_SplineFFT(5, 64, T),
        dt,
    ) < 1e-6
    T = Double64
    dt = T(1)/11
#    @time @test test_translation((100, 120), [Lagrange(15, T), Lagrange(15, T)], 11) < 1e-14
    @time @test test_translation((100, 120), Lagrange(15, T), Lagrange(15, T), dt) < 1e-14
    @time @test test_translation(
        (128, 64),
        B_SplineLU(15, 128, T),
        B_SplineLU(15, 64, T),
        dt,
    ) < 1e-18
    # @time @test test_translation(
    #     (128, 64),
    #     [B_SplineLU(15, 128, T), B_SplineLU(15, 64, T)],
    #     11,
    # ) < 1e-17
    @time @test test_translation(
        (128, 64),
        B_SplineFFT(15, 128, T),
        B_SplineFFT(15, 64, T),
        dt,
    ) < 1e-17
    # @time @test test_translation(
    #     (128, 64),
    #     [B_SplineFFT(15, 128, T), B_SplineFFT(15, 64, T)],
    #     11,
    # ) < 1e-17
    T = BigFloat
    dt = T(1)/11
#    @time @test test_translation((100, 90), [Lagrange(21, T), Lagrange(21, T)], 11) < 1e-17
@time @test test_translation((100, 90), Lagrange(25, T), Lagrange(25, T), dt) < 1e-20
@time @test test_translation(
        (128, 64),
        B_SplineLU(25, 128, T),
        B_SplineLU(25, 64, T),
        dt,
    ) < 1e-22
    # @time @test test_translation(
    #     (128, 64),
    #     [B_SplineLU(21, 128, T), B_SplineLU(21, 64, T)],
    #     11,
    # ) < 1e-22
    @time @test test_translation(
        (128, 64),
        B_SplineFFT(25, 128, T),
        B_SplineFFT(25, 64, T),
        dt,
    ) < 1e-22
   # @time @test test_translation(
    #     (128, 64),
    #     [B_SplineFFT(21, 128, T), B_SplineFFT(21, 64, T)],
    #     11,
    # ) < 1e-22

    # T = Float64
    # @time @test test_translation(
    #     (32, 38, 30),
    #     [B_SplineLU(9, 32, T), Lagrange(9, T), Lagrange(9, T)],
    #     11,
    # ) < 1e-4
    # @time @test test_translation(
    #     (38, 32, 30),
    #     [Lagrange(9, T), B_SplineLU(9, 32, T), Lagrange(9, T)],
    #     11,
    # ) < 1e-4
    # @time @test test_translation(
    #     (30, 38, 32),
    #     [Lagrange(9, T), Lagrange(9, T), B_SplineLU(9, 32, T)],
    #     11,
    # ) < 1e-4
    # @time @test test_translation(
    #     (32, 64, 30),
    #     [B_SplineLU(9, 32, T), B_SplineFFT(9, 64, T), Lagrange(9, T)],
    #     11,
    # ) < 1e-4
    # @time @test test_translation(
    #     (30, 32, 64),
    #     [Lagrange(9, T), B_SplineLU(9, 32, T), B_SplineFFT(9, 64, T)],
    #     11,
    # ) < 1e-4
    # @time @test test_translation(
    #     (32, 30, 64),
    #     [B_SplineLU(9, 32, T), Lagrange(9, T), B_SplineFFT(9, 64, T)],
    #     11,
    # ) < 1e-4

end
