using SemiLagrangian:
    advection!,
    split_advection!,
    AdvectionOperator,
    strangsplit,
    interpolate!,
    LocalLagrange,
    LagrangeFFT,
    CubicSpline,
    Spectral,
    Lagrange,
    BSplineLU,
    UniformMesh

struct TestAdvectionInterpolant end

function SemiLagrangian.interpolate!(output, input, shift, ::TestAdvectionInterpolant)
    return output .= input .+ shift
end

include(joinpath(@__DIR__, "..", "examples", "landau-damping-split.jl"))

@testset "Landau operator example" begin
    nx, nv, dt, nt = 16, 24, 0.1, 3
    meshx = UniformMesh(0.0, 4pi, nx)
    meshv = UniformMesh(-6.0, 6.0, nv)
    x, v = meshx.points, meshv.points
    f = [exp(-vi^2 / 2) * (1 + 0.001cos(0.5xi)) / sqrt(2pi) for xi in x, vi in v]
    interp_x, interp_v = Spectral(nx), Spectral(nv)
    field = compute_e(meshx, compute_rho(meshv, f))
    expected_energy = [sum(abs2, field) * step(meshx)]
    for _ in 1:nt
        advection!(f, meshx, v, dt / 2, interp_x; dim = 1)
        field = compute_e(meshx, compute_rho(meshv, f))
        push!(expected_energy, sum(abs2, field) * step(meshx))
        advection!(f, meshv, field, dt, interp_v; dim = 2)
        advection!(f, meshx, v, dt / 2, interp_x; dim = 1)
    end
    time, energy = landau_split(nx, nv, dt, nt)
    @test time ≈ [0.0, 0.05, 0.15, 0.25]
    @test energy ≈ expected_energy
end

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

    for interpolant in (LocalLagrange(7), LagrangeFFT(n, 7), CubicSpline(n), Spectral(n))
        shifted = similar(initial)
        interpolate!(shifted, initial, 0, interpolant)
        @test shifted ≈ initial
        interpolate!(shifted, initial, 0.25f0, interpolant)
        float32_shifted = copy(shifted)
        interpolate!(shifted, initial, 0.25, interpolant)
        @test shifted ≈ float32_shifted atol = 1.0e-7
        @test_throws ArgumentError interpolate!(shifted, initial, missing, interpolant)
        @test_throws MethodError interpolate!(shifted, interpolant, initial, shift)

        interpolate!(shifted, initial, shift, interpolant)
        result = copy(initial)
        @test advection!(result, mesh, T(0.35), T(1) / n, interpolant) === result
        @test result ≈ shifted
    end

    @testset "self-consistent Landau splitting" begin
        nx, nv = 16, 24
        meshx = UniformMesh(0.0, 4pi, nx)
        meshv = UniformMesh(-6.0, 6.0, nv)
        x, v = meshx.points, meshv.points
        landau_initial = [exp(-vi^2 / 2) * (1 + 0.01cos(0.5xi)) / sqrt(2pi) for xi in x, vi in v]
        dt = 0.1
        coefficients = strangsplit(dt)

        for (interp_x, interp_v) in (
                (Spectral(nx), Spectral(nv)),
                (LagrangeFFT(nx, 7), LagrangeFFT(nv, 7)),
                (LocalLagrange(7), LocalLagrange(7)),
                (CubicSpline(nx), CubicSpline(nv)),
            )
            reference = copy(landau_initial)
            actual = copy(landau_initial)
            transposed = similar(actual, nv, nx)
            position = AdvectionOperator(meshx, v, interp_x; dim = 2)
            initial_field = compute_e(meshx, compute_rho(meshv, actual))

            for _ in 1:3
                advection!(reference, meshx, v, dt / 2, interp_x; dim = 1)
                reference_field = compute_e(meshx, compute_rho(meshv, reference))
                advection!(reference, meshv, reference_field, dt, interp_v; dim = 2)
                advection!(reference, meshx, v, dt / 2, interp_x; dim = 1)

                advection!(actual, meshx, v, coefficients[1], interp_x; dim = 1)
                field = compute_e(meshx, compute_rho(meshv, actual))
                @test field ≈ reference_field
                @test maximum(abs.(field .- initial_field)) > 1.0e-8
                permutedims!(transposed, actual, (2, 1))
                velocity = AdvectionOperator(meshv, field, interp_v; dim = 1)
                split_advection!(transposed, velocity, position, @view(coefficients[2:3]))
                permutedims!(actual, transposed, (2, 1))
                @test actual ≈ reference
            end
        end
    end

    custom_result = copy(initial)
    @test advection!(custom_result, mesh, T(0.35), T(1) / n, TestAdvectionInterpolant()) ===
        custom_result
    @test custom_result ≈ initial .+ shift

    spline_interp = BSplineLU(4, n)
    velocities = [0.2 * sin(2pi * x) for x in mesh.points]
    variable_expected = similar(initial)
    for i in eachindex(initial)
        line = similar(initial)
        interpolate!(line, initial, -velocities[i] * T(0.1) / step(mesh), spline_interp)
        variable_expected[i] = line[i]
    end
    variable_result = copy(initial)
    @test advection!(variable_result, mesh, velocities, T(0.1), spline_interp) ===
        variable_result
    @test variable_result ≈ variable_expected

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
