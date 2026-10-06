using LinearAlgebra
using SemiLagrangian:
    Polynomial,
    derivative,
    Hermite,
    HermiteCache,
    hermite_interp,
    hermite_interp!,
    PrecalHermite,
    L,
    Lprim,
    K,
    H,
    bplus,
    bminus,
    _getlagrangecoefficients

lagrange_polynomial_reference(k, order, origin) =
    Polynomial(_getlagrangecoefficients(k, order, origin))

@testset "Hermite one-shot interpolation" begin
    x = [0.0, 0.15, 0.4, 0.8, 1.1, 1.7, 2.0]
    f(x) = x^3 - 2x^2 + 4x - 1
    df(x) = 3x^2 - 4x + 4
    y = f.(x)
    dy = df.(x)
    targets = reshape([0.0, 0.25, 0.7, 1.3, 2.0], 1, :)
    expected = f.(targets)

    @test hermite_interp(x, y, dy, targets) ≈ expected
    @test hermite_interp(x, y, dy, 0.7) ≈ f(0.7)

    output = similar(expected)
    @test hermite_interp!(output, x, y, dy, targets) === output
    @test output ≈ expected

    cache = HermiteCache(x)
    @test hermite_interp(cache, y, dy, targets) ≈ expected
    cached_output = similar(expected)
    @test hermite_interp!(cached_output, cache, y, dy, targets) === cached_output
    @test cached_output ≈ expected

    y2 = sin.(x)
    dy2 = cos.(x)
    @test hermite_interp(cache, y2, dy2, targets) ≈ sin.(targets) atol = 5e-4

    @test_throws DomainError hermite_interp(cache, y, dy, -0.1)
    @test_throws DimensionMismatch hermite_interp(cache, y[1:end-1], dy, targets)
    @test_throws DimensionMismatch hermite_interp!(zeros(2), cache, y, dy, targets)
    @test_throws ArgumentError HermiteCache([0.0, 1.0, 1.0])
end

function getbp(i, rp, sp)
    res = big(1 // 1)
    for j = rp:sp
        if j != i
            res /= (i - j)
            if j != 0
                res *= (-j)
            end
        end
    end
    return res
end

function test_precalhermite(ord)
    ph = PrecalHermite(ord)
    d = div(ord, 2)
    @test ph.rplus == -d
    @test ph.splus == d + 1
    @test ph.rminus == -d - 1
    @test ph.sminus == d

    sbp = 0 // 1
    sbm = 0 // 1

    for i = (-d):(d+1)
        basis = lagrange_polynomial_reference(i + d, ord, -d)
        @test L(ph, i) == basis
        @test Lprim(ph, i) == derivative(L(ph, i))(i)
        @test K(ph, i) == basis^2 * Polynomial([-i, 1 // 1])
        @test H(ph, i) ==
              basis^2 * (1 - 2 * Lprim(ph, i) * Polynomial([-i, 1 // 1]))
        if i != 0
            @test bplus(ph, i) == getbp(i, ph.rplus, ph.splus)
            @test bminus(ph, -i) == -getbp(i, ph.rplus, ph.splus)
        end
        sbp += bplus(ph, i)
        sbm += bminus(ph, -i)
    end

    @test sbp == 0
    @test sbm == 0

    if ord == 3
        @test ph.bplus == [-1 // 3, -1 // 2, 1, -1 // 6]
    elseif ord == 5
        @test ph.bplus == [1 // 20, -1 // 2, -1 // 3, 1 // 1, -1 // 4, 1 // 30]
    end
end
function test_base_hermite(order)
    herm = Hermite(order, Rational{BigInt})
    ord = div(order, 2) + 1
    tab = rationalize.(BigInt, rand(ord + 1), tol = 1 / 1000000)
    fct = Polynomial(tab)
    dec = div(order, 2)
    for i = 1:3
        value = rationalize.(BigInt, rand(), tol = 1 / 1000000)
        res = 0 // 1
        for j = 0:order
            res += herm.tabfct[j+1](value) * fct(j - dec)
        end
        @test res == fct(value)
    end
end
function test_base_hermite2d(order)
    herm = Hermite(order, Rational{BigInt})
    ord = div(order, 2) + 1
    tab = rationalize.(BigInt, rand(ord + 1, ord + 1), tol = 1 / 1000000)
    fct = Pol2(tab)
    dec = div(order, 2)
    for i = 1:1
        val_x = rationalize.(BigInt, rand(), tol = 1 / 1000000)
        val_y = rationalize.(BigInt, rand(), tol = 1 / 1000000)
        res = 0 // 1
        for j = 0:order, k = 0:order
            res += herm.tabfct[j+1](val_x) * herm.tabfct[k+1](val_y) * fct(j - dec, k - dec)
        end
        resf = fct(val_x, val_y)
        #        println("order = $order res-fct=$(res-resf)")
        @test res == resf
    end
end
@time @testset "Hermite precal" begin
    for ord = 3:2:21
        test_precalhermite(ord)
    end
end
@time @testset "test base hermite interpolation" begin
    for order = 5:4:41
        test_base_hermite(order)
    end
end
@time @testset "test base hermite interpolation2d" begin
    for order = 5:4:41
        test_base_hermite2d(order)
    end
end
