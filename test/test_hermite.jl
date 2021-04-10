using LinearAlgebra
using Polynomials
using SemiLagrangian: Hermite



function test_base_hermite(order)
    herm = Hermite(order, Rational{BigInt})
    tab = rationalize.(BigInt, rand(order), tol = 1 / 1000000)
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
    tab = rationalize.(BigInt, rand(order, order), tol = 1 / 1000000)
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
@time @testset "test base hermite interpolation" begin
    for ord = 5:4:41
        test_base_hermite(ord)
    end
end
@time @testset "test base hermite interpolation2d" begin
    for ord = 5:4:25
        test_base_hermite2d(ord)
    end
end
