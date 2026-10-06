using Test

@testitem "Spline interpolation" begin
    nx = 100
    alpha = 0.2
    u = Float64[cos(2π * (i - 1) / nx) for i = 1:nx]
    u_out = zeros(nx)
    expected = Float64[cos(2π * ((i - 1) + alpha) / nx) for i = 1:nx]

    for order in [2, 4, 6, 8]

        interpolant = PeriodicBSpline(nx, order)
        interpolate!(u_out, interpolant, u, 0.0)
        @test maximum(abs.(u_out - u)) < 1e-14

        interpolate!(u_out, interpolant, u, alpha)
        tol = max(10.0 / 10^(2order), 1e-14)
        @test maximum(abs.(u_out - expected)) < tol

    end
end

@testitem "B-Splines basis" begin
    p = 3
    biatx = SemiLagrangian.uniform_bsplines_eval_basis(p, 0.0)
    @test biatx[1] ≈ 1/6
    @test biatx[2] ≈ 2/3
    @test biatx[3] ≈ 1/6
end

@testitem "Arbitrary degree B-spline interpolation" begin
    using SemiLagrangian: SplineCache, spline_interp, spline_interp!

    x = [0.0, 0.1, 0.28, 0.46, 0.7, 0.91, 1.0]
    queries = reshape([0.0, 0.05, 0.22, 0.51, 0.83, 1.0], 2, 3)
    for order in 1:5
        f(t) = t^order
        values = f.(x)
        expected = f.(queries)
        cache = SplineCache(x; order)

        @test spline_interp(cache, values, queries) ≈ expected atol = 1e-10
        output = similar(expected)
        @test spline_interp!(output, cache, values, queries) === output
        @test output ≈ expected atol = 1e-10
        @test spline_interp(x, values, queries; order) ≈ expected atol = 1e-10
        @test spline_interp(cache, values, 0.37) ≈ f(0.37) atol = 1e-10
    end

    cache = SplineCache(x; order = 3)
    @test spline_interp(cache, cos.(x), queries) ≈
          spline_interp(x, cos.(x), queries; order = 3)
    @test_throws DomainError spline_interp(cache, x, -0.1)
    @test_throws DimensionMismatch spline_interp(cache, x[1:end-1], queries)
    @test_throws ArgumentError SplineCache(x; order = length(x))
end
