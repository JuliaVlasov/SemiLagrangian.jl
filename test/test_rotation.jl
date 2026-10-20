using DoubleFloats
using LinearAlgebra
using SemiLagrangian:
    advection!,
    magicsplit,
    UniformMesh,
    AbstractInterpolation,
    Lagrange,
    Hermite,
    BSplineLU,
    BSplineFFT

function exact_rotation!(f, mesh_x, mesh_y, t)
    s, c = sincos(t)
    for (i, x) in enumerate(mesh_x.points), (j, y) in enumerate(mesh_y.points)
        x0, y0 = c * x - s * y, s * x + c * y
        f[i, j] = exp(-2 * (x0^2 + (y0 + typeof(t)(6 // 5))^2))
    end
    return f
end

function test_rotation(
    sz,
    interp_x::AbstractInterpolation{T},
    interp_y::AbstractInterpolation{T},
    nsteps,
) where {T}
    mesh_x = UniformMesh(T(-5), T(5), sz[1])
    mesh_y = UniformMesh(T(-6), T(4.5), sz[2])
    dt = T(2big(pi) / nsteps)
    data = zeros(T, sz)
    exact = similar(data)
    exact_rotation!(data, mesh_x, mesh_y, zero(T))

    for n in 1:nsteps
        for (i, coefficient) in enumerate(magicsplit(dt))
            if isodd(i)
                advection!(data, mesh_x, mesh_y.points, coefficient, interp_x; dim = 1)
            else
                advection!(data, mesh_y, -mesh_x.points, coefficient, interp_y; dim = 2)
            end
        end
        exact_rotation!(exact, mesh_x, mesh_y, dt * n)
    end
    return norm(data .- exact, Inf)
end

@testset "directional rotation" begin
    T = Float64
    @test test_rotation((2000, 1022), Lagrange(5, T), Lagrange(5, T), 11) < 1e-3
    @test test_rotation((2000, 1022), Hermite(5, T), Hermite(5, T), 11) < 1e-3
    @test test_rotation(
        (128, 256),
        BSplineLU(5, 128, T),
        BSplineLU(5, 256, T),
        11,
    ) < 1e-3
    @test test_rotation(
        (128, 256),
        BSplineFFT(5, 128, T),
        BSplineFFT(5, 256, T),
        11,
    ) < 1e-3

    T = Double64
    @test test_rotation((200, 220), Lagrange(15, T), Lagrange(15, T), 11) < 1e-7
    @test test_rotation(
        (128, 256),
        BSplineLU(15, 128, T),
        BSplineLU(15, 256, T),
        11,
    ) < 1e-8

    T = BigFloat
    @test test_rotation((200, 220), Lagrange(25, T), Lagrange(25, T), 11) < 1e-9
    @test test_rotation(
        (128, 256),
        BSplineFFT(25, 128, T),
        BSplineFFT(25, 256, T),
        11,
    ) < 1e-9
end
