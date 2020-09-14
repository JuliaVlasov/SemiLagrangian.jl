include("../src/advection.jl")
include("../src/spline.jl")
include("../src/bspline.jl")
include("../src/bsplinelu.jl")
include("../src/bsplinefft.jl")
include("../src/lagrange.jl")

using LinearAlgebra
using Test

function test_interp(interp::InterpolationType{Rational{BigInt}, iscirc}, sz) where {iscirc}

    @time @testset "test interpolation  $(get_type(interp))" begin    
        fct = if iscirc
            getbspline(big"3",0)
        else
            x -> x^3-x^2-x//6+1//4
        end
        dec=(big"10113"//big"10240")*1
        mesh=collect(4big.(0:(sz-1)))//sz
        deb = fct.(mesh)
 #       println("deb=$deb")
        ref=fct.(mesh .+ 4dec//sz)
        res = zeros(Rational{BigInt}, sz)
        interpolate!(missing, res, deb, dec, interp)
        println("norm=$(float(norm(res-ref,Inf)))")
        diff = res-ref
        # for i=1:sz
        #     println("i=$i diff=$(convert(Float64,diff[i])), $(diff[i])")
        # end
        @test res == ref
    end
end

lag=Lagrange(Rational{BigInt},3,iscirc=false)
test_interp(lag,128)
bsp=B_SplineLU(3,128,Rational{BigInt}; iscirc=true)
test_interp(bsp,128)
# bsp=B_SplineLU(3,128,Rational{BigInt}; iscirc=false)
# test_interp(bsp,128)
