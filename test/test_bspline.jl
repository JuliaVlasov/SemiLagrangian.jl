
include("../src/advection.jl")
include("../src/spline.jl")
include("../src/bspline.jl")
include("../src/bsplinelu.jl")
include("../src/bsplinefft.jl")

using LinearAlgebra
using Test

function test_spline(order, prec)
    setprecision(prec) do
        @time @testset "test spline object of order $order" begin
           p0 = getbspline(order,0)
            p3 = getbspline(order,3)
            @test p3 == decal(p0,3)
            @test isapprox(sum(p0.(range(big"1.0",length=order+1))),1, atol=1e-60)
            @test sum(p0.(range(one(Rational{BigInt}),length=order+1))) == 1             
            for i = 0:order
                x = rand(BigFloat)
                res0=p0(x+i)
                res3=p3(x+i+3)
                @test isapprox( res0, res3,atol=1e-60)
                @test order == Polynomials.degree(p0[i])
                pol1 = p0[i]
                pol2 = p0[i+1]
                for j=1:order
                    @test pol1(i+1) == pol2(i+1)
                    pol1 = derivative(pol1)
                    pol2 = derivative(pol2)
                end
            end
        end
    end
end
for order = 1:15
    test_spline( order, order < 20 ? 256 : 512 )
end
function nb_diff_tol(f1, f2, tol)
    nb=0
    for i=1:size(f1,1)
        if abs(f1[i]-f2[i])> tol
  #          println("indice=$i diff=$(abs(f1[i]-f2[i]))")
            nb+=1
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
function bsplinerec(p, j, x::T) where{T}

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
        tab_x=[1821//10000, 1//1234, 7677//8999, big"123456789"//big"234567890", 0//1]
        for ord=1:9, x in tab_x
            @test bsplinerec(3,0,x) == getbspline(3,0)(x)
        end
    end

end
       
function test_interpolation(type::DataType, order, iscirc::Bool, n, nb,  tol, islu::Bool)
    
    sp = if (islu)
            B_SplineLU(order, n, zero(type); iscirc=iscirc)
    else
        B_SplineFFT(order, n, zero(type))
    end
    coef = convert(type, iscirc ? 1 : big"1.111")
#    fct(v,n) = exp( -cos(2big(pi)*coef*v/n)^2)
#    fct(v,n) = exp(-(75*(v-n/2)/n)^2)
 #   fct(v,n) = exp( -(cos(2big(pi)*coef*v/n))^2)
    fct(v,n) = cos(2big(pi)*coef*v/n)
 fp = fct.(convert.(type,(collect(1:n))),n)
    fi = zeros(type, n)
    value = convert(type,
#    big"0.000000000000000000000000000000000000001") 
    -big"0.385713901112334905767655546588878787878787887874441132001118762519")
#    big"0.385713901112334905767655546588878787878787887874441619187615524132001118762519")
    for i=1:nb
        fi .= fp
#        fpref = fct.(convert.(type,(collect(1:n))) .+ i*value,n)
        fpref = fct.((1:n) .+ i*value,n)
        interpolate!(fp, fi, value, sp)
        # for i=1:number
        #     println("i=$i norm=$(norm(fpref[i]-fp[i]))")
        # end
 #       if i == n
        if i%100 == 0
            println("i=$i norm=$(norm(fpref-fp)), nb=$(nb_diff_tol(fpref,fp,1e-40))")
        end
        println("i=$i norm=$(norm(fpref-fp))")
 #       end
        @test isapprox(fpref, fp, atol=tol)
        # println("spline=$(sp.bspline)")
        # break
    end
end

function test_interpolation_2d(type::DataType, order, iscirc::Bool, n,  tol)
    
    sp = B_SplineLU(order, n, zero(type); iscirc=iscirc)
    coef = convert(type, iscirc ? 1 : big"1.111")
#    fct(v,n) = exp( -cos(2big(pi)*coef*v/n)^2)
    fct(x,y,n) = cos((coef*2big(pi)/n)*(x+y))
    fp = [fct(convert.(type,x), convert.(type, y), n) for x=1:n, y=1:n]
    fi = zeros(type, n)
    buf = zeros(type, n)
    value_x = convert(type, 
#    big"0.0385713901112334905767655546588878787878787887874441619187615524132001118762519")
    big"0.385713901112334905767655546588878787878787887874441619187615524132001118762519")
    value_y = convert(type, 
    #    big"0.0385713901112334905767655546588878787878787887874441619187615524132001118762519")
        big"0.13901112334905767655546588878787878787887874441619187615524132001118762519")
    k=0
    fpref = [fct.(x+k*value_x, y+k*value_y, n) for x=1:n, y=1:n]
    println("k=$k norm=$(norm(fpref-fp))")
    @test isapprox(fpref, fp, atol=tol)

    for k=1:n, l=1:2
        for i=1:n
            fi .= fp[:,i]
            interpolate!( buf, fi, (value_x, value_y)[l], sp)
            fp[:,i] .= buf
        end
        fp=transpose(fp)
        if l == 2
            fpref = [fct.(x+k*value_x, y+k*value_y, n) for x=1:n, y=1:n]
            # for i=1:n
            #     println("k=$k i=$i norm=$(norm(fpref[:,i]-fp[:,i]))")
            # end
            println("k=$k norm=$(norm(fpref-fp))")
            @test isapprox(fpref, fp, atol=tol)
        end
    end

end

test_bspline()

@testset "test interpolation bspline" begin
    # test_interpolation(BigFloat, 11, true, 100, 1e-7)
    # test_interpolation(BigFloat, 21, true, 100, 1e-15)
    # test_interpolation(BigFloat, 41, true, 100, 1e-20)
    # @time test_interpolation(BigFloat, 21, true, 2^14, 100, 1e-10, false)
    # @time test_interpolation(BigFloat, 21, true, 2^14, 100, 1e-10, true)
    #@time test_interpolation(BigFloat, 9, true, 2^8, 100, 1, false)
    @time test_interpolation(BigFloat, 9, true, 2^8, 100, 1, true)
    # test_interpolation_2d(BigFloat, 27, true, 100, 1e-20)
end


