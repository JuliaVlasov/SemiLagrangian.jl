using DoubleFloats
using LinearAlgebra

using SemiLagrangian:
    Advection1d,
    sizeall,
    Advection1dData,
    getdata,
    advection!,
    UniformMesh,
    gettranslationvar1d,
    AbstractInterpolation,
    interpolate!,
    Lagrange,
    B_SplineLU,
    B_SplineFFT
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

function test_translation1d(
    sz::NTuple{2,Int},
    interp_sp::AbstractInterpolation{T},
    interp_v::AbstractInterpolation{T},
    nbdt::Int,
) where {T}
@time begin
    spmin, spmax, nsp = T(-5), T(5), sz[1]
    vmin, vmax, nv = -T(5.55), T(5), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)
end
println("trace11d")
@time begin
    dt = T(1)
    adv = Advection1d((mesh_sp,), (mesh_v,), (interp_sp,), (interp_v,), dt)
    tabref = zeros(T, sz)

    v1 = T(big"0.83545655467782872872782870029282982828737872878776717190927267611111")
    v2 = T(-big"0.945678101929276765616767176761671771766717828781828998101092981877817176")
end
println("trace21d")
@time begin 
    exact!(tabref, (v1, v2), T(0))

    pvar = gettranslationvar1d((v1, v2))

    advdata = Advection1dData(adv, tabref, pvar)

    diffmax = 0
    data = getdata(advdata)
end
println("trace31d")
    for ind = 1:nbdt
        while advection!(advdata)
        end
        #        GC.gc()
        exact!(tabref, (v1, v2), T(ind))
        diff = norm(data .- tabref, Inf)
        #        @show v1, v2, ind, diff
        diffmax = max(diffmax, diff)
    end
    println(
        "test_translation1d sz=$sz interp=$interp_sp, $interp_v nbdt=$nbdt diffmax=$diffmax",
    )
    return diffmax

end
function test_translation1d(
    sz::NTuple{N,Int},
    interp_t::Vector{I},
    nbdt::Int,
) where {T,N, I<:AbstractInterpolation{T}}

    dt = T(1)

    vall = ntuple(x -> T(10rand(BigFloat) - 5), N)

    vall =  T(1)/10 .* vall

    @show vall

    tabref = zeros(T, sz)
    exact!(tabref, vall, T(0))

    data = copy(tabref)

    buf = zeros(T, sz)

    # dec1 = fill(v1, sz)
    # dec2 = fill(v2, sz)


    diffmax = 0
    for ind = 1:nbdt
        interpolate!(buf, data, (i, _) ->vall[i], interp_t)
        #        interpolate!(buf, data, (dec1, dec2), interp)
        copyto!(data, buf)
        exact!(tabref, vall, T(ind))
        diff = norm(data .- tabref, Inf)
        diffmax = max(diffmax, diff)
        @show ind, diff
    end
    println("test_translation1d sz=$sz interp=$interp_t nbdt=$nbdt diffmax=$diffmax")
    return diffmax

end


@testset "test translation1d" begin
    T = Float64
    # @time @test test_translation1d((150, 220), [Lagrange(9, T), Lagrange(9, T)], 11) < 1e-6
    # @time @test test_translation1d((50, 30, 40), map(x -> Lagrange(5+2x, T),1:3), 11) < 1e-4
    # @time @test test_translation1d((100, 120), [Lagrange(5, T), Lagrange(5, T)], 11) < 1e-5
    # @time @test test_translation1d((128, 64), [B_SplineLU(9, 128, T), B_SplineLU(9, 64, T)], 11) < 1e-6
    @time @test test_translation1d(
        (200, 300),
        Lagrange(5, T),
        Lagrange(5, T),
        11,
    ) < 1e-3
    @time @test test_translation1d(
        (128, 64),
        B_SplineLU(5, 128, T),
        B_SplineLU(5, 64, T),
        11,
    ) < 1e-6
    @time @test test_translation1d(
        (128, 64),
        B_SplineFFT(5, 128, T),
        B_SplineFFT(5, 64, T),
        11,
    ) < 1e-6
    T = Double64
#    @time @test test_translation1d((100, 120), [Lagrange(15, T), Lagrange(15, T)], 11) < 1e-14
    @time @test test_translation1d((100, 120), Lagrange(15, T), Lagrange(15, T), 11) < 1e-14
    @time @test test_translation1d(
        (128, 64),
        B_SplineLU(15, 128, T),
        B_SplineLU(15, 64, T),
        11,
    ) < 1e-18
    # @time @test test_translation1d(
    #     (128, 64),
    #     [B_SplineLU(15, 128, T), B_SplineLU(15, 64, T)],
    #     11,
    # ) < 1e-17
    @time @test test_translation1d(
        (128, 64),
        B_SplineFFT(15, 128, T),
        B_SplineFFT(15, 64, T),
        11,
    ) < 1e-17
    # @time @test test_translation1d(
    #     (128, 64),
    #     [B_SplineFFT(15, 128, T), B_SplineFFT(15, 64, T)],
    #     11,
    # ) < 1e-17
    T = BigFloat
#    @time @test test_translation1d((100, 90), [Lagrange(21, T), Lagrange(21, T)], 11) < 1e-17
    @time @test test_translation1d((100, 90), Lagrange(25, T), Lagrange(25, T), 11) < 1e-20
    @time @test test_translation1d(
        (128, 64),
        B_SplineLU(25, 128, T),
        B_SplineLU(25, 64, T),
        11,
    ) < 1e-22
    # @time @test test_translation1d(
    #     (128, 64),
    #     [B_SplineLU(21, 128, T), B_SplineLU(21, 64, T)],
    #     11,
    # ) < 1e-22
    @time @test test_translation1d(
        (128, 64),
        B_SplineFFT(25, 128, T),
        B_SplineFFT(25, 64, T),
        11,
    ) < 1e-22
    # @time @test test_translation1d(
    #     (128, 64),
    #     [B_SplineFFT(21, 128, T), B_SplineFFT(21, 64, T)],
    #     11,
    # ) < 1e-22

    # T = Float64
    # @time @test test_translation1d((32, 38, 30),[ B_SplineLU(9, 32, T), Lagrange(9, T), Lagrange(9,T)], 11) < 1e-4
    # @time @test test_translation1d((38, 32, 30),[ Lagrange(9,T), B_SplineLU(9, 32, T), Lagrange(9, T)], 11) < 1e-4
    # @time @test test_translation1d((30, 38, 32),[ Lagrange(9, T), Lagrange(9,T), B_SplineLU(9, 32, T)], 11) < 1e-4
    # @time @test test_translation1d((32, 64, 30),[ B_SplineLU(9, 32, T), B_SplineFFT(9, 64, T), Lagrange(9,T)], 11) < 1e-4
    # @time @test test_translation1d((30, 32, 64),[ Lagrange(9,T), B_SplineLU(9, 32, T), B_SplineFFT(9, 64, T), ], 11) < 1e-4
    # @time @test test_translation1d((32, 30, 64),[ B_SplineLU(9, 32, T), Lagrange(9,T), B_SplineFFT(9, 64, T), ], 11) < 1e-4

end
