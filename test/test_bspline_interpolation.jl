using Test

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
