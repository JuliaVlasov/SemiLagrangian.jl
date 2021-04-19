using LinearAlgebra
using Polynomials
using SemiLagrangian: Hermite, PrecalHermite, L, Lprim, K, H, bplus, bminus, _getpolylagrange


function getbp(i,rp,sp)
    res=big(1//1)
    for j=rp:sp
        if j != i
            res /= (i-j)
            if j != 0
                res *= (-j)
            end
        end
    end
    return res
end

function test_precalhermite(ord)
    ph = PrecalHermite(ord)
    d = div(ord,2)
    @test ph.rplus == -d
    @test ph.splus == d+1
    @test ph.rminus == -d-1
    @test ph.sminus == d

    sbp=0//1
    sbm=0//1

    for i=-d:d+1
        @test L(ph,i) == _getpolylagrange(i+d,ord,-d)
        @test Lprim(ph, i) == derivative(L(ph,i))(i)
        @test K(ph, i) == _getpolylagrange(i+d,ord,-d)^2 *Polynomial([-i, 1//1])
        @test H(ph, i) == _getpolylagrange(i+d,ord,-d)^2 * (1-2*Lprim(ph,i)*Polynomial([-i,1//1]))
        if i != 0
            @test bplus(ph, i) == getbp(i,ph.rplus,ph.splus)
            @test bminus(ph, -i) == -getbp(i, ph.rplus, ph.splus)
        end
        sbp += bplus(ph,i)
        sbm += bminus(ph, -i)
    end

    @test sbp == 0
    @test sbm == 0

    if ord == 3
        @test ph.bplus == [-1//3, -1//2, 1, -1//6]
    elseif ord == 5
        @test ph.bplus == [1//20,-1//2,-1//3,1//1,-1//4,1//30]
    end

end
function test_base_hermite(order)
    herm = Hermite(order, Rational{BigInt})
    ord = div(order,2)+1
    tab = rationalize.(BigInt, rand(ord+1), tol = 1 / 1000000)
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
    ord = div(order,2)+1
    tab = rationalize.(BigInt, rand(ord+1, ord+1), tol = 1 / 1000000)
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
    for ord=3:2:21
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

