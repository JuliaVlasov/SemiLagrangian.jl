using LinearAlgebra
using SemiLagrangian:
    UniformMesh,
    AbstractInterpolation,
    Lagrange,
    interpolate!,
    gettabmod,
    CachePrecal

function test_swirling(
    sz::NTuple{2,Int},
    interp::AbstractVector{I},
    nsteps::Int,
) where {T,I<:AbstractInterpolation{T}}
    mesh_x = UniformMesh(T(0), T(1), sz[1])
    mesh_y = UniformMesh(T(0), T(1), sz[2])
    dt = T(1.5) / nsteps
    velocity_x = zeros(T, sz)
    velocity_y = zeros(T, sz)
    data = zeros(T, sz)

    for i in 1:sz[1], j in 1:sz[2]
        x, y = mesh_x.points[i], mesh_y.points[j]
        velocity_x[i, j] = -cos(T(pi) * x)^2 * sin(2T(pi) * y) / step(mesh_x)
        velocity_y[i, j] = cos(T(pi) * y)^2 * sin(2T(pi) * x) / step(mesh_y)
        data[i, j] = ((1 - x)^2 + (1 - y)^2 < 0.8) ? 1 : 0
    end

    tabmod = gettabmod.(sz)
    cache = CachePrecal(interp, one(T))
    buffer = similar(data)
    for n in 1:nsteps
        coefficient = dt * cos(T(pi) * (n - 1) / nsteps)
        interpolate!(
            buffer,
            data,
            x -> (coefficient * velocity_x[x], coefficient * velocity_y[x]),
            interp;
            tabmod = tabmod,
            cache = cache,
        )
        copyto!(data, buffer)
    end
    return data
end

@testset "variable characteristic interpolation" begin
    result = test_swirling((100, 100), [Lagrange(9), Lagrange(9)], 50)
    @test all(isfinite, result)
end
