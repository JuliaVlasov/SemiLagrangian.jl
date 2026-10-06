using FastInterpolations

@testset "FastInterpolations periodic cubic interpolation" begin
    n = 64
    xi = 1.0:n
    alpha = 0.25
    u = sin.(2π .* (0:n-1) ./ n)
    expected = sin.(2π .* ((0:n-1) .+ alpha) ./ n)

    result = cubic_interp(
        xi,
        u,
        xi .+ alpha;
        bc=PeriodicBC(endpoint=:exclusive),
    )

    @test result ≈ expected atol=1e-5
    @test cubic_interp(xi, u, xi; bc=PeriodicBC(endpoint=:exclusive)) ≈ u atol=1e-12
end
