using LinearAlgebra
using SemiLagrangian:
    Polynomial,
    derivative,
    Hermite,
    ExtrapolationBC,
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

@testset "Hermite boundary conditions" begin
    extrapolation_interp = Hermite(5; bc = ExtrapolationBC())
    @test extrapolation_interp isa Hermite{Float64, ExtrapolationBC, 5}
end

function getbp(i, rp, sp)
    res = big(1 // 1)
    for j in rp:sp
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

    for i in (-d):(d + 1)
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

    return if ord == 3
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
    for i in 1:3
        value = rationalize.(BigInt, rand(), tol = 1 / 1000000)
        res = 0 // 1
        for j in 0:order
            res += herm.tabfct[j + 1](value) * fct(j - dec)
        end
        @test res == fct(value)
    end
    return
end
function test_base_hermite2d(order)
    herm = Hermite(order, Rational{BigInt})
    ord = div(order, 2) + 1
    tab = rationalize.(BigInt, rand(ord + 1, ord + 1), tol = 1 / 1000000)
    fct = Pol2(tab)
    dec = div(order, 2)
    for i in 1:1
        val_x = rationalize.(BigInt, rand(), tol = 1 / 1000000)
        val_y = rationalize.(BigInt, rand(), tol = 1 / 1000000)
        res = 0 // 1
        for j in 0:order, k in 0:order
            res += herm.tabfct[j + 1](val_x) * herm.tabfct[k + 1](val_y) * fct(j - dec, k - dec)
        end
        resf = fct(val_x, val_y)
        #        println("order = $order res-fct=$(res-resf)")
        @test res == resf
    end
    return
end
@time @testset "Hermite precal" begin
    for ord in 3:2:21
        test_precalhermite(ord)
    end
end
@time @testset "test base hermite interpolation" begin
    for order in 5:4:41
        test_base_hermite(order)
    end
end
@time @testset "test base hermite interpolation2d" begin
    for order in 5:4:41
        test_base_hermite2d(order)
    end
end
