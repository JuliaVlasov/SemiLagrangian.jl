using Test
using SemiLagrangian

include(joinpath(@__DIR__, "..", "examples", "landau-damping-poisson.jl"))

@testset "Poisson advection" begin
    for nx in (15, 16)
        nv = 24
        meshx = UniformMesh(0.0, 4pi, nx)
        meshv = UniformMesh(-6.0, 6.0, nv)
        operator = PoissonAdvectionOperator(meshx, meshv)
        density = 1 .+ 0.1 .* cos.(0.5 .* meshx.points)
        f = repeat(density ./ meshv.width, 1, nv)
        original = copy(f)
        field = update_field!(operator, f)
        @test f == original
        @test operator.rho ≈ density .- 1 atol = 1.0e-14
        @test field ≈ 0.2 .* sin.(0.5 .* meshx.points) atol = 1.0e-14
        @test abs(sum(field)) < 1.0e-14
        fill!(f, 1.0)
        @test maximum(abs, update_field!(operator, f)) < 1.0e-14
        if iseven(nx)
            for i in 1:nx
                f[i, :] .= isodd(i) ? 1.0 : -1.0
            end
            @test maximum(abs, update_field!(operator, f)) < 1.0e-14
        end
        @test_throws DimensionMismatch update_field!(operator, zeros(nx, nv + 1))
    end

    nx, nv, dt = 16, 24, 0.1
    meshx = UniformMesh(0.0, 4pi, nx)
    meshv = UniformMesh(-6.0, 6.0, nv)
    initial = [
        (1 + 0.001cos(0.5x)) * exp(-v^2 / 2) / sqrt(2pi)
            for x in meshx.points, v in meshv.points
    ]
    for interpolant in (Spectral(nv), LocalLagrange(7), BSplineLU(4, nv))
        operator = PoissonAdvectionOperator(meshx, meshv, interpolant)
        actual = copy(initial)
        reference = copy(initial)
        update_field!(operator, reference)
        advection!(reference, meshv, operator.field, dt, interpolant; dim = 2)
        @test advection!(actual, operator, dt) === actual
        @test actual ≈ reference
    end

    operator = PoissonAdvectionOperator(meshx, meshv)
    position = AdvectionOperator(meshx, meshv.points, Spectral(nx); dim = 1)
    actual, reference = copy(initial), copy(initial)
    coefficients = strangsplit(dt)
    expected_energy = zeros(4)
    update_field!(operator, initial)
    expected_energy[1] = sum(abs2, operator.field) * step(meshx)
    for it in 1:3
        advection!(reference, meshx, meshv.points, dt / 2, position.interpolant; dim = 1)
        update_field!(operator, reference)
        expected_energy[it + 1] = sum(abs2, operator.field) * step(meshx)
        advection!(reference, meshv, operator.field, dt, Spectral(nv); dim = 2)
        advection!(reference, meshx, meshv.points, dt / 2, position.interpolant; dim = 1)
        split_advection!(actual, position, operator, coefficients)
        @test actual ≈ reference
    end
    time, energy = landau_poisson(nx, nv, dt, 3)
    @test time ≈ [0, 0.05, 0.15, 0.25]
    @test energy ≈ expected_energy

    # Measure after compilation, including a complete Poisson and velocity stage.
    advection!(actual, operator, dt)
    allocated_stage(data, op, dt) = @allocated advection!(data, op, dt)
    allocated_stage(actual, operator, dt)
    @test allocated_stage(actual, operator, dt) == 0
    reverse_reference = copy(initial)
    reverse_actual = copy(initial)
    advection!(reverse_reference, operator, dt)
    advection!(reverse_reference, position, dt)
    split_advection!(reverse_actual, operator, position, standardsplit(dt))
    @test reverse_actual ≈ reverse_reference
    @test_throws ArgumentError split_advection!(actual, position, operator, Float64[])
    @test_throws DimensionMismatch PoissonAdvectionOperator(meshx, meshv, Spectral(nv + 1))
    @test_throws ArgumentError PoissonAdvectionOperator(UniformMesh(1.0, 0.0, nx), meshv)
end
