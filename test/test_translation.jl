using DoubleFloats
using LinearAlgebra
using SemiLagrangian:
    advection!,
    strangsplit,
    standardsplit,
    triplejumpsplit,
    UniformMesh,
    Lagrange,
    BSplineLU,
    BSplineFFT,
    Hermite

function exact_translation!(f, velocities, t)
    T = typeof(t)
    for ind in CartesianIndices(f)
        f[ind] = exp(
            -sum(
                sin(2T(pi) * (ind.I[i] - velocities[i] * t) / size(f, i) +
                    T(big"0.56") * i^2) for i in eachindex(velocities)
            ),
        )
    end
    return f
end

function test_translation(sz, interp_sp, interp_v, dt; tc = strangsplit(dt))
    T = typeof(dt)
    mesh_sp = UniformMesh(T(-5), T(5), sz[1])
    mesh_v = UniformMesh(T(-5.55), T(5), sz[2])
    velocities = (
        T(big"0.83545655467782872872782870029282982828737872878776717190927267611111"),
        T(big"0.945678101929276765616767176761671771766717828781828998101092981877817176"),
    )
    exact = zeros(T, sz)
    exact_translation!(exact, velocities, zero(T))
    data = copy(exact)
    velocity_sp = fill(velocities[1] * step(mesh_sp), sz[2])
    velocity_v = fill(velocities[2] * step(mesh_v), sz[1])

    for n in 1:round(Int, one(T) / dt)
        for (i, coefficient) in enumerate(tc)
            if isodd(i)
                advection!(data, mesh_sp, velocity_sp, coefficient, interp_sp; dim = 1)
            else
                advection!(data, mesh_v, velocity_v, coefficient, interp_v; dim = 2)
            end
        end
        exact_translation!(exact, velocities, T(n) * dt)
    end
    return norm(data .- exact, Inf)
end

@testset "directional translation" begin
    T = Float64
    dt = T(1) / 11
    for split in (strangsplit(dt), triplejumpsplit(dt), standardsplit(dt))
        @test test_translation((200, 300), Lagrange(5, T), Lagrange(5, T), dt; tc = split) <
              1e-3
        @test test_translation((200, 300), Hermite(5, T), Hermite(5, T), dt; tc = split) <
              1e-3
    end
    @test test_translation(
        (128, 64),
        BSplineLU(5, 128, T),
        BSplineLU(5, 64, T),
        dt,
    ) < 1e-6
    @test test_translation(
        (128, 64),
        BSplineFFT(5, 128, T),
        BSplineFFT(5, 64, T),
        dt,
    ) < 1e-6

    T = Double64
    dt = T(1) / 11
    @test test_translation((100, 120), Lagrange(15, T), Lagrange(15, T), dt) < 1e-14
    @test test_translation(
        (128, 64),
        BSplineLU(15, 128, T),
        BSplineLU(15, 64, T),
        dt,
    ) < 1e-18

    T = BigFloat
    dt = T(1) / 11
    @test test_translation((100, 90), Lagrange(25, T), Lagrange(25, T), dt) < 1e-20
    @test test_translation(
        (128, 64),
        BSplineFFT(25, 128, T),
        BSplineFFT(25, 64, T),
        dt,
    ) < 1e-22
end
