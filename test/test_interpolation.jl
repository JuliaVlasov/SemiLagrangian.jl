include("../src/advection.jl")
include("../src/spline.jl")
include("../src/bspline.jl")
include("../src/bsplinelu.jl")
include("../src/bsplinefft.jl")
include("../src/lagrange.jl")
include("../src/interpolation.jl")

using LinearAlgebra
using Test

function test_interp(interp::InterpolationType{Rational{BigInt}, iscirc}, sz) where {iscirc}

    @time @testset "test interpolation  $(get_type(interp))" begin    
        fct = if iscirc
            order=3
            getbspline(big(order),0)
        else
            order=0
            x -> x^3-x^2-x//6+1//4
        end
        dec=(big"1"//big"10240")*1
        mesh=collect((order+1)*big.(0:(sz-1)))//sz
        deb = fct.(mesh)
 #       println("deb=$deb")
        ref=fct.(mesh .+ (order+1)*dec//sz)
        res = zeros(Rational{BigInt}, sz)
        interpolate!(res, deb, dec, interp)
        println("norm=$(float(norm(res-ref,Inf)))")
        diff = res-ref
        # for i=1:sz
        #     println("i=$i diff=$(convert(Float64,diff[i])), $(diff[i])")
        # end
        fl = !iscirc || isbspline(interp)
        println("fl=$fl isbspline=$(isbspline(interp))")
        if fl
            @test res == ref
        else
            @test size(filter( x->x==0,res-ref),1)==120
        end
    end
end

lag=Lagrange(Rational{BigInt},3,iscirc=false)
test_interp(lag,128)
lag=Lagrange(Rational{BigInt},3,iscirc=true)
test_interp(lag,128)
bsp=B_SplineLU(3,128,Rational{BigInt}; iscirc=true)
test_interp(bsp,128)
# bsp=B_SplineLU(3,128,Rational{BigInt}; iscirc=false)
# test_interp(bsp,128)
