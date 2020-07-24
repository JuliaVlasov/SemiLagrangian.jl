using Test
using LinearAlgebra
include("../src/lagrange2d.jl")

function test_baselagrange2d(order)
    N = (order > 10) ? BigInt : Int64
    origin=-div(order,2)
    @testset "Lagrange interpolation 2d order=$order" begin
        for i=0:order, j=0:order
            pl = _getpolylagrange2d(i, j, order, origin, N)
            if order == 13
                println("i=$i j=$j type pl = $(typeof(pl))")
            end
            for x=origin:origin+order, y=origin:origin+order
                @test pl(N(x),N(y)) == ((x==i+origin && y==j+origin) ? 1 : 0)
            end
        end
    end
end

function test_lagrange2d(type::DataType, order, iscirc::Bool, number,  tol)
    @testset "Lagrange interpolation 2d function type=$type order=$order iscirc=$iscirc" begin
        lag = Lagrange2d(type, order; iscirc=iscirc)
#        println("lag=$lag")
        n = number
        # whennot circuar coef != 1 set the function unperiodic
        coef = iscirc ? big"1." : big"1.111"
        fct(x,y,n) = cos(2big(pi)*coef*x/n+big"0.43")*sin(2big(pi)*coef*y/n+big"0.2")
        tf = [fct(i,j,n) for i=1:n, j=1:n]
        value_x = big"0.456231986"
        value_y = big"0.3676956102"
        # value_x = big"0.01111561987109"
        # value_y = big"0.011217009811"
        normax=0.0
        for i=1:number, j=1:number
            pol = polinterpol(lag, tf, (i,j))
            ref=fct(i+value_x,j+value_y,n)
            res=pol(value_x, value_y)
            normax = max(normax, norm(ref-res))
#            println("i=$i norm=$(norm(res-ref,Inf))")
            # println("i=$i j=$j norm=$(norm(res-ref,Inf)) type(pol)=$(typeof(pol)) pol=$pol")
            # println("fct(i,j,n)=$(fct(i,j,n))")
            # @test isapprox(ref, res, atol=tol)
        end
        println("normax=$normax")
    end
end

function test_interpolation(type::DataType, order, iscirc::Bool, number, granularity,  tol, nb=1)
    
    lag = Lagrange2d(type, order; iscirc=iscirc, granularity=granularity)
    coef = iscirc ? 1. : 1.111
    n = number
    fct(x,y,n) = exp( -cos(2big(pi)*coef*v/n)^2+ sin(2big(pi)*coef*y/n+big"0.2"))
    fp = [fct(i,j,n) for i=1:n, j=1:n]
    fi = zeros(type, number, number)
    value_x = big"0.485713901"
    value_y = big"0.185713901"
    for i=1:nb
        fi .= fp
        fpref = [fct(x+i*value_x, y+i*value_y,n) for x=1:n, y=1:n]
        interpolate!(missing, fp, fi, (value_x, value_y), lag)
        println("i=$i granularity=$granularity norm = $(norm(fpref-fp,Inf))")
        # for i=1:number
        #     println("i=$i norm=$(norm(fpref[i]-fp[i]))")
        # end
        @test isapprox(fpref, fp, atol=tol)
    end
end
for i=3:7
    test_baselagrange2d(i)
end
# test_baselagrange2d(17)

test_lagrange2d(BigFloat, 5,true, 100, 1e-10)
test_interpolation(BigFoat, 5, true, 100, 1, 1e-7, 1000)
# test_lagrange2d(BigFloat, 5,true, 100, 1e-50)
# test_lagrange(BigFloat, 17,false, 1000, 1e-38)
# test_lagrange(BigFloat, 23,false, 1000, 1e-50)
# @testset "Lagrange interpolation order=17 iscirc=true" begin
#     for i= 1:25
#         @time test_interpolation(BigFloat, 17,true, 1000, i, 1e-20)
#     end
# end
# @time test_interpolation(BigFloat, 17,true, 200, 1, 1e-10, 200 )
# test_interpolation(BigFloat, 23,true, 1000, 12, 1e-40)
# test_interpolation(BigFloat, 17,false, 1000, 10, 1e-38)
# test_interpolation(BigFloat, 23,false, 1000, 12, 1e-40)






