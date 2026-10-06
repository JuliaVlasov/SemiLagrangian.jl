using LinearAlgebra
using SemiLagrangian:
    Polynomial,
    Lagrange,
    _getlagrangecoefficients,
    _c,
    LagrangeCache,
    lagrange_interp,
    lagrange_interp!

struct Pol2{T}
    tab::Array{T,2}
end
function (fct::Pol2{T})(x, y) where {T}
    res = 0
    pc_x = [T(x)^(i - 1) for i = 1:size(fct.tab, 1)]
    pc_y = [T(y)^(i - 1) for i = 1:size(fct.tab, 2)]
    for i = 1:size(fct.tab, 1), j = 1:size(fct.tab, 2)
        res += fct.tab[i, j] * pc_x[i] * pc_y[j]
    end
    return res
end

function test_base_lagrange(order)
    lag = Lagrange(order, Rational{BigInt})
    tab = rationalize.(BigInt, rand(order + 1), tol = 1 / 1000000)
    fct = Polynomial(tab)
    dec = div(order, 2)
    for i = 1:3
        value = rationalize.(BigInt, rand(), tol = 1 / 1000000)
        res = 0 // 1
        for j = 0:order
            res += lag.tabfct[j+1](value) * fct(j - dec)
        end
        @test res == fct(value)
    end
end
function test_base_lagrange2d(order)
    lag = Lagrange(order, Rational{BigInt})
    tab = rationalize.(BigInt, rand(order + 1, order + 1), tol = 1 / 1000000)
    fct = Pol2(tab)
    dec = div(order, 2)
    for i = 1:1
        val_x = rationalize.(BigInt, rand(), tol = 1 / 1000000)
        val_y = rationalize.(BigInt, rand(), tol = 1 / 1000000)
        res = 0 // 1
        for j = 0:order, k = 0:order
            res += lag.tabfct[j+1](val_x) * lag.tabfct[k+1](val_y) * fct(j - dec, k - dec)
        end
        resf = fct(val_x, val_y)
        #        println("order = $order res-fct=$(res-resf)")
        @test res == resf
    end
end
@testset "Lagrange basis coefficient construction" begin
    for order = 0:15, origin in (-div(order, 2), 0), k = 0:order
        coefficients = _getlagrangecoefficients(k, order, origin)
        for x in (big(-1 // 3), big(0 // 1), big(2 // 5))
            reference = prod(
                ((x - (l + origin)) / (k - l) for l = 0:order if l != k);
                init = one(x),
            )
            @test evalpoly(x, coefficients) == reference
        end
    end

    for order = 0:12, k = 0:order
        coefficients = _getlagrangecoefficients(k, order, 0)
        integral = sum(
            coefficients[m] * ((isodd(m) ? 1 : -1) // m) for m in eachindex(coefficients)
        )
        @test _c(k, order) == integral
    end
end

@testset "Lagrange one-shot interpolation" begin
    x = [0.0, 0.15, 0.4, 0.8, 1.1, 1.7, 2.0]
    y = x .^ 3 .- 2x .^ 2 .+ 4x .- 1
    targets = reshape([0.0, 0.25, 0.7, 1.3, 2.0], 1, :)
    expected = targets .^ 3 .- 2targets .^ 2 .+ 4targets .- 1

    result = lagrange_interp(x, y, targets; order = 3)
    @test result ≈ expected

    output = similar(result)
    @test lagrange_interp!(output, x, y, targets; order = 3) === output
    @test output ≈ expected
    @test lagrange_interp(x, y, 0.7; order = 3) ≈ 0.7^3 - 2 * 0.7^2 + 4 * 0.7 - 1
    @test lagrange_interp([0, 2, 3], [0, 4, 9], [1]; order = 1) ≈ [2.0]

    cache = LagrangeCache(x; order = 3)
    @test lagrange_interp(cache, y, targets) ≈ expected
    cached_output = similar(output)
    @test lagrange_interp!(cached_output, cache, y, targets) === cached_output
    @test cached_output ≈ expected
    y2 = cos.(x)
    @test lagrange_interp(cache, y2, targets) ≈ lagrange_interp(x, y2, targets; order = 3)
    @test_throws DimensionMismatch lagrange_interp(cache, y[1:end-1], targets)

    @test_throws DomainError lagrange_interp(x, y, -0.1; order = 3)
    @test_throws DimensionMismatch lagrange_interp!(zeros(2), x, y, targets; order = 3)
    @test_throws ArgumentError lagrange_interp([0.0, 1.0, 1.0], [1.0, 2.0, 3.0], [0.5]; order = 1)
    @test_throws ArgumentError lagrange_interp(x, y, [0.5]; order = -1)
end

@time @testset "test base interpolation" begin
    for ord = 3:27
        test_base_lagrange(ord)
    end
end
@time @testset "test base interpolation2d" begin
    for ord = 3:21
        test_base_lagrange2d(ord)
    end
end
