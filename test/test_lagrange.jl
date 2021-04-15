using LinearAlgebra
using Polynomials
using SemiLagrangian: Lagrange

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
    tab = rationalize.(BigInt, rand(order), tol = 1 / 1000000)
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
    tab = rationalize.(BigInt, rand(order, order), tol = 1 / 1000000)
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
