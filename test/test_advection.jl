using SemiLagrangian: advection!, interpolate!, Lagrange, UniformMesh

@testset "directional advection" begin
    T = Float64
    n = 64
    mesh = UniformMesh(T(0), T(1), n)
    interp = Lagrange(5, T)
    initial = sin.(2T(pi) .* mesh.points)
    expected = similar(initial)
    shift = -T(0.35)

    interpolate!(expected, initial, shift, interp)
    result = copy(initial)
    @test advection!(result, mesh, T(0.35), T(1) / n, interp) === result
    @test result ≈ expected

    mesh_y = UniformMesh(T(-1), T(1), 32)
    data = [sin(2T(pi) * x) * cos(T(pi) * y) for x in mesh.points, y in mesh_y.points]
    original = copy(data)
    velocity = collect(range(T(-0.5), T(0.5); length = size(data, 1)))
    out = similar(data)
    for i in axes(data, 1)
        interpolate!(
            view(out, i, :),
            view(original, i, :),
            -velocity[i] * T(0.1) / step(mesh_y),
            interp,
        )
    end
    @test advection!(data, mesh_y, velocity, T(0.1), interp; dim = 2) === data
    @test data ≈ out

    pointwise_velocity = [T(0.2) + T(0.1) * x for x in mesh.points, _ in mesh_y.points]
    pointwise_result = copy(original)
    @test advection!(
        pointwise_result,
        mesh,
        pointwise_velocity,
        T(0.05),
        interp;
        dim = 1,
    ) === pointwise_result
    @test all(isfinite, pointwise_result)

    @test_throws DimensionMismatch advection!(copy(initial), mesh_y, one(T), T(0.1), interp)
    @test_throws DimensionMismatch advection!(copy(initial), mesh, ones(3), T(0.1), interp)
    @test_throws ArgumentError advection!(copy(initial), mesh, one(T), T(0.1), interp; dim = 2)
end
