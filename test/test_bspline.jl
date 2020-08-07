include("../src/bspline.jl")
include("../src/matspline.jl")

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

function test_interpolation(type::DataType, order, iscirc::Bool, n,  tol)
    
    sp = BSplineNew(order, n, zero(type); iscirc=iscirc)
    coef = convert(type, iscirc ? 1 : 1.111)
#    fct(v,n) = exp( -cos(2big(pi)*coef*v/n)^2)
    fct(v,n) = cos((coef*2big(pi)/n)*v)
    fp = fct.(convert.(type,(collect(1:n))),n)
    fi = zeros(type, n)
    value = convert(type, 
#    big"0.0385713901112334905767655546588878787878787887874441619187615524132001118762519")
    -big"1.385713901112334905767655546588878787878787887874441619187615524132001118762519")
    for i=1:n
        fi .= fp
#        fpref = fct.(convert.(type,(collect(1:n))) .+ i*value,n)
        fpref = fct.((1:n) .+ i*value,n)
        interpolate!(missing, fp, fi, value, sp)
        # for i=1:number
        #     println("i=$i norm=$(norm(fpref[i]-fp[i]))")
        # end
        if i == n
            println("i=$i norm=$(norm(fpref-fp))")
        end
        @test isapprox(fpref, fp, atol=tol)
    end
end

function test_interpolation_2d(type::DataType, order, iscirc::Bool, n,  tol)
    
    sp = BSplineNew(order, n, zero(type); iscirc=iscirc)
    coef = convert(type, iscirc ? 1 : 1.111)
#    fct(v,n) = exp( -cos(2big(pi)*coef*v/n)^2)
    fct(x,y,n) = cos((coef*2big(pi)/n)*(x+y))
    fp = [fct(convert.(type,x), convert.(type, y), n) for x=1:n, y=1:n]
    fi = zeros(type, n)
    buf = zeros(type, n)
    value_x = convert(type, 
#    big"0.0385713901112334905767655546588878787878787887874441619187615524132001118762519")
    -big"1.385713901112334905767655546588878787878787887874441619187615524132001118762519")
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
            interpolate!(missing, buf, fi, (value_x, value_y)[l], sp)
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

@testset "test interpolation bspline" begin
    test_interpolation(BigFloat, 11, true, 100, 1e-2)
    test_interpolation(BigFloat, 21, true, 100, 1e-3)
    test_interpolation(BigFloat, 41, true, 100, 1e-2)
    test_interpolation_2d(BigFloat, 27, true, 100, 1e-2)
end
