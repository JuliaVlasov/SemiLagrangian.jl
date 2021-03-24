

using LinearAlgebra
using Polynomials

using SemiLagrangian:
    getbspline, B_SplineLU, get_kl_ku, B_SplineFFT, interpolate!, decal, get_precal

function test_spline(order, prec)
    setprecision(prec) do
        @time @testset "test spline object of order $order" begin
            p0 = getbspline(order, 0)
            p3 = getbspline(order, 3)
            @test p3 == decal(p0, 3)
            @test isapprox(sum(p0.(range(big"1.0", length = order + 1))), 1, atol = 1e-60)
            @test sum(p0.(range(one(Rational{BigInt}), length = order + 1))) == 1
            for i = 0:order
                x = rand(BigFloat)
                res0 = p0(x + i)
                res3 = p3(x + i + 3)
                @test isapprox(res0, res3, atol = 1e-60)
                @test order == Polynomials.degree(p0[i])
                pol1 = p0[i]
                pol2 = p0[i+1]
                for j = 1:order
                    @test pol1(i + 1) == pol2(i + 1)
                    pol1 = derivative(pol1)
                    pol2 = derivative(pol2)
                end
            end
        end
    end
end
for order = 1:15
    test_spline(order, order < 20 ? 256 : 512)
end
function nb_diff_tol(f1, f2, tol)
    nb = 0
    for i = 1:size(f1, 1)
        if abs(f1[i] - f2[i]) > tol
            #          println("indice=$i diff=$(abs(f1[i]-f2[i]))")
            nb += 1
        end
    end
    return nb
end
"""
    bsplinerec(p, j, x)

Return the value at x in [0,1[ of the B-spline with
integer nodes of degree p with support starting at j.
Implemented recursively using the de Boor's recursion formula
using the [De Boor's Algorithm](https://en.wikipedia.org/wiki/De_Boor%27s_algorithm)

Note that this function has an exponential complexity depending to the order, so it is only there for testing.

```math
B_{i,0}(x) := \\left\\{
\\begin{matrix}
1 & \\mathrm{if}  \\quad t_i ≤ x < t_{i+1} \\\\
0 & \\mathrm{otherwise}
\\end{matrix}
\\right.
```

```math
B_{i,p}(x) := \\frac{x - t_i}{t_{i+p} - t_i} B_{i,p-1}(x)
+ \\frac{t_{i+p+1} - x}{t_{i+p+1} - t_{i+1}} B_{i+1,p-1}(x).
```

"""
function bsplinerec(p, j, x::T) where {T}

    if p == 0
        if j == 0
            return one(T)
        else
            return zero(T)
        end
    else
        w = (x - j) / p
        w1 = (x - j - 1) / p
    end
    (w * bsplinerec(p - 1, j, x) + (1 - w1) * bsplinerec(p - 1, j + 1, x))

end

function test_bspline()

    @time @testset "test verify bspline" begin
        tab_x = [
            1821 // 10000,
            1 // 1234,
            7677 // 8999,
            big"123456789" // big"234567890",
            0 // 1,
        ]
        for ord = 1:10, x in tab_x
            if ord >8
                p1 = Polynomial([big(1//7),2//7,3//7])
                p2 = Polynomial([-big(6//7), 2//7, 3//7])
            else
                p1 = Polynomial([1//7,2//7,3//7])
                p2 = Polynomial([-6//7, 2//7, 3//7])
            end
            for i = 0:(ord-1)
                bsp = getbspline(ord, i)
                b1 = bsp*p1
                b2 = bsp*p2
                bspbis = b1 - b2
                v = bsplinerec(ord, i, x)
                @test  v == bsp(x)
                @test v == bspbis(x)
            end
        end
    end

end



@testset "test get_kl_ku" begin
    kl, ku = get_kl_ku(5)
    @test kl == 2 && ku == 2

    kl, ku = get_kl_ku(6)
    @test kl == 2 && ku == 3
    #    @test kl == 3 && ku == 2
end

test_bspline()
