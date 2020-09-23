include("../src/advection.jl")
include("../src/spline.jl")
include("../src/bspline.jl")
include("../src/bsplinelu.jl")
include("../src/bsplinefft.jl")
include("../src/lagrange.jl")
include("../src/interpolation.jl")

using LinearAlgebra
using Test

function test_interp(interp::InterpolationType{Rational{BigInt}, iscirc}, dec,  sz) where {iscirc}

    @time @testset "test interpolation  $(get_type(interp)) dec=$dec" begin    
        fct = if iscirc
            order=3
            getbspline(big(order),0)
        else
            order=0
            x -> x^3-x^2-x//6+1//4
        end
        mesh=collect((order+1)*big.(0:(sz-1)))//sz
        deb = fct.(mesh)
 #       println("deb=$deb")
        ref=if iscirc
            decf = mod(dec,1//1)
            decint = div(dec,1//1)
            circshift(fct.(mesh .+ (order+1)*decf//sz), (-decint,))
        else
            fct.(mesh .+ (order+1)*dec//sz)
        end
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
function test_interp(interp, sz)
    for i=0:5
        println("test_interp i=$i")
        test_interp(interp, big"1235"//10240 +i, sz)
    end
end
function test_interpfloat(interp::InterpolationType{T, iscirc}, dec::T, sz, tol, fct) where {T,iscirc}

    @time @testset "test interpolation  $(get_type(interp)) order=$(get_order(interp))" begin    
        
        mesh=T.(collect(big.(0:(sz-1)))/sz)
        deb = fct.(mesh)
 #       println("deb=$deb")
        ref=fct.(mesh .+ dec/sz)
        res = zeros(T, sz)
        interpolate!(res, deb, dec, interp)
        println("norm=$(float(norm(res-ref,Inf)))")
        diff = res-ref
        # for i=1:sz
        #     println("i=$i diff=$(convert(Float64,diff[i])), $(diff[i])")
        # end
        @test isapprox(ref, res, atol=tol)
    end
end

lag=Lagrange(Rational{BigInt},3,iscirc=false)
test_interp(lag,big"3"//1024,128)
lag=Lagrange(Rational{BigInt},3,iscirc=true)
test_interp(lag,128)
bsp=B_SplineLU(3,128,Rational{BigInt}; iscirc=true)
test_interp(bsp,128)
# bsp = B_SplineLU(3,128,Float64; iscirc=true)
# fct(x) = cos(2pi*x+0.25)
# test_interpfloat(bsp,128,1e-6,fct)
# bsp = B_SplineLU(31,128,BigFloat; iscirc=true)
fct(x) = cos(2big(pi)*x+big"0.25")
# test_interpfloat(bsp,128,1e-50, fct)
# fct(x) = exp(-260*(x-0.4)^2)
# test_interpfloat(bsp,128,1e-18, fct)
# bsp=B_SplineLU(3,128,Rational{BigInt}; iscirc=false)
# test_interp(bsp,128)
lag=Lagrange(BigFloat,3,iscirc=true)
test_interpfloat(lag,big"0.34111111111191919191",128,1e-5, fct)
test_interpfloat(lag,-big"0.24111111111191919191",128,1e-5, fct)
test_interpfloat(lag,big"5.34111111111191919191",128,1e-5, fct)
test_interpfloat(lag,-big"3.34111111111191919191",128,1e-5, fct)
lag=Lagrange(BigFloat,21,iscirc=true)
test_interpfloat(lag,big"0.34111111111191919191",128,1e-20, fct)
test_interpfloat(lag,-big"0.24111111111191919191",128,1e-20, fct)
test_interpfloat(lag,big"5.34111111111191919191",128,1e-20, fct)
test_interpfloat(lag,-big"3.34111111111191919191",128,1e-20, fct)
lag=Lagrange(BigFloat,4,iscirc=true)
test_interpfloat(lag,big"0.34111111111191919191",128,1e-5, fct)
test_interpfloat(lag,-big"0.24111111111191919191",128,1e-5, fct)
test_interpfloat(lag,big"5.34111111111191919191",128,1e-5, fct)
test_interpfloat(lag,-big"3.34111111111191919191",128,1e-5, fct)
lag=Lagrange(BigFloat,22,iscirc=true)
test_interpfloat(lag,big"0.34111111111191919191",128,1e-20, fct)
test_interpfloat(lag,-big"0.24111111111191919191",128,1e-20, fct)
test_interpfloat(lag,big"5.34111111111191919191",128,1e-20, fct)
test_interpfloat(lag,-big"3.34111111111191919191",128,1e-20, fct)

bsp=B_SplineLU(3,128,BigFloat; iscirc=true)
test_interpfloat(bsp,big"0.34111111111191919191",128,1e-7, fct)
test_interpfloat(bsp,-big"0.24111111111191919191",128,1e-7, fct)
test_interpfloat(bsp,big"5.34111111111191919191",128,1e-7, fct)
test_interpfloat(bsp,-big"3.34111111111191919191",128,1e-7, fct)
bsp=B_SplineLU(21,128,BigFloat; iscirc=true)
test_interpfloat(bsp,big"0.34111111111191919191",128,1e-40, fct)
test_interpfloat(bsp,-big"0.24111111111191919191",128,1e-40, fct)
test_interpfloat(bsp,big"5.34111111111191919191",128,1e-40, fct)
test_interpfloat(bsp,-big"3.34111111111191919191",128,1e-40, fct)
bsp=B_SplineLU(4,129,BigFloat; iscirc=true)
test_interpfloat(bsp,big"0.34111111111191919191",129,1e-7, fct)
test_interpfloat(bsp,-big"0.24111111111191919191",129,1e-7, fct)
test_interpfloat(bsp,big"5.34111111111191919191",129,1e-7, fct)
test_interpfloat(bsp,-big"3.34111111111191919191",129,1e-7, fct)
bsp=B_SplineLU(22,129,BigFloat; iscirc=true)
test_interpfloat(bsp,big"0.34111111111191919191",129,1e-40, fct)
test_interpfloat(bsp,-big"0.24111111111191919191",129,1e-40, fct)
test_interpfloat(bsp,big"5.34111111111191919191",129,1e-40, fct)
test_interpfloat(bsp,-big"3.34111111111191919191",129,1e-40, fct)
fct(x)=exp(-260*(x-1.4)^2)
bsp=B_SplineLU(3,128,BigFloat; iscirc=true)
test_interpfloat(bsp,big"0.34111111111191919191",128,1e-4, fct)
test_interpfloat(bsp,-big"0.24111111111191919191",128,1e-4, fct)
test_interpfloat(bsp,big"5.34111111111191919191",128,1e-4, fct)
test_interpfloat(bsp,-big"3.34111111111191919191",128,1e-4, fct)
bsp=B_SplineLU(21,128,BigFloat; iscirc=true)
test_interpfloat(bsp,big"0.34111111111191919191",128,1e-18, fct)
test_interpfloat(bsp,-big"0.24111111111191919191",128,1e-18, fct)
test_interpfloat(bsp,big"5.34111111111191919191",128,1e-15, fct)
test_interpfloat(bsp,-big"3.34111111111191919191",128,1e-15, fct)
bsp=B_SplineLU(4,129,BigFloat; iscirc=true)
test_interpfloat(bsp,big"0.34111111111191919191",129,1e-4, fct)
test_interpfloat(bsp,-big"0.24111111111191919191",129,1e-4, fct)
test_interpfloat(bsp,big"5.34111111111191919191",129,1e-4, fct)
test_interpfloat(bsp,-big"3.34111111111191919191",129,1e-4, fct)
bsp=B_SplineLU(22,129,BigFloat; iscirc=true)
test_interpfloat(bsp,big"0.34111111111191919191",129,1e-18, fct)
test_interpfloat(bsp,-big"0.24111111111191919191",129,1e-18, fct)
test_interpfloat(bsp,big"5.34111111111191919191",129,1e-15, fct)
test_interpfloat(bsp,-big"3.34111111111191919191",129,1e-15, fct)
