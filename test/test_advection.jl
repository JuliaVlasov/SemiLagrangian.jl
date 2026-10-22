using SemiLagrangian:
    advection!,
    split_advection!,
    AdvectionOperator,
    strangsplit,
    interpolate!,
    Lagrange,
    UniformMesh

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

    dt = T(0.1)
    velocity_x = collect(range(-T(0.5), T(0.5); length = size(data, 2)))
    operator_x = AdvectionOperator(mesh, velocity_x, interp; dim = 1)
    velocity_y = collect(range(-T(0.5), T(0.5); length = n))
    operator_y = AdvectionOperator(mesh_y, velocity_y, interp; dim = 2)
    coefficients = strangsplit(dt)
    expected_split = copy(original)
    for (stage, coefficient) in enumerate(coefficients)
        if isodd(stage)
            advection!(expected_split, mesh, velocity_x, coefficient, interp; dim = 1)
        else
            advection!(
                expected_split,
                mesh_y,
                velocity_y,
                coefficient,
                interp;
                dim = 2,
            )
        end
    end
    split_result = copy(original)
    @test split_advection!(split_result, operator_x, operator_y, coefficients) ===
          split_result
    @test split_result ≈ expected_split
    @test_throws ArgumentError split_advection!(
        copy(original),
        operator_x,
        operator_y,
        T[],
    )
end
