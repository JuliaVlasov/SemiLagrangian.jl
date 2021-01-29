
include("../src/mesh.jl")
include("../src/advection.jl")
include("../src/lagrange.jl")
# include("../src/lagrange2d.jl")
# include("../src/advection2d.jl")
include("../src/bspline.jl")
include("../src/bsplinelu.jl")
include("../src/bsplinefft.jl")
include("../src/interpolation.jl")

# using Images
using LinearAlgebra
using Test
"""

   exact(tf, mesh1, mesh2)

   Julia function to compute exact solution "

```math
\frac{d f}{dt} +  (y \frac{df}{delta1} - x \frac{df}{delta2}) = 0 
```

"""
function exact(tf::T, mesh1::UniformMesh{T}, mesh2::UniformMesh{T}) where {T}

    f = zeros(T,(mesh1.length,mesh2.length))
    for (i, x) in enumerate(mesh1.points), (j, y) in enumerate(mesh2.points)
        xn = cos(tf) * x - sin(tf) * y
        yn = sin(tf) * x + cos(tf) * y
#       f[i,j] = (xn-0.3)^2+(yn+0.5)^2 < 0.03 ? 1.0 : 0.0
# f[i,j] = exp(-5*(cos(xn-big"0.65")^2+sin(yn+big"0.5")^2))
        f[i,j] = exp(-13*((xn)^2+(yn+big"1.2")^2))
#       f[i,j] = exp(-12*((xn)^2+(yn+big"1.8")^2))
#       f[i,j] = exp(-(xn-1)*(xn-1)*10)*exp(-(yn)*(yn)*10)
#       f[i,j] = exp(-(sin(xn)+0.4)^2)*exp(-(cos(yn)-0.5)^2)
#       f[i,j] = exp(-sin(xn-0.3)*sin(xn-0.3)/0.6)*exp(-cos(yn+0.5)*cos(yn+0.5)/0.6)
    end

    f

end

" Function to compute error "
function error1(f, f_exact)
    i_max=0
    j_max=0
    v_max=0
    i_min=0
    j_min=0
    v_min=1.0
    for i=1:size(f,1), j=1:size(f,2)
        v = abs(f[i,j]-f_exact[i,j])
        if v > v_max
            i_max, j_max, v_max = i, j, v
        end
        if v < v_min
            i_min, j_min, v_min = i, j, v
        end
    end
#    println("i_max=$i_max j_max=$j_max v_max=$v_max")
#    println("i_min=$i_min j_min=$j_min v_min=$v_min")
    return v_max
    # maximum(abs.(f .- f_exact))
end

# function savefile( f, par, str)
#     img=zeros(RGB,size(f))
#     k = -log.(f)
#     k .-= minimum(k)
#     maxk = maximum(k)
#     k ./= maxk

#     borne = min(0.5, par/maxk)

#     for i=1:size(f,1), j=1:size(f,2)
#         v = k[i,j]
#         b = v > borne ? (v-borne)/(1-borne) : 0
#         r = 1 -v
#         g = v < borne ? (borne-v)/borne : 0
#         img[i,j] = RGB(r, g, b)
#     end
#     Images.save(str, img)
# end



function rotation_2d(
    tf::T, 
    nt, 
    mesh1::UniformMesh{T}, 
    mesh2::UniformMesh{T}, 
    interp::InterpolationType
) where {T}

    dt = tf/nt

    # n1 = mesh1.length
    # x1min, x1max = mesh1.start, mesh1.stop
    delta1 = mesh1.step

    # n2 = mesh2.length
    # x2min, x2max = mesh2.start, mesh2.stop
    delta2 = mesh2.step

#    println("delta1=$delta1 delta2=$delta2")

    f  = zeros(T,(mesh1.length,mesh2.length))
    f .= exact(zero(T), mesh1, mesh2)
    ft = zeros(T,(mesh1.length,mesh2.length))
#    transpose!(ft, f)

#    println(f)

    
    adv_x1 = Advection( mesh1, interp )
    adv_x2 = Advection( mesh2, interp )

    v1 = collect(mesh2.points)
    v2 =  - collect(mesh1.points)


    # println("v1=$v1")
    # println("v2=$v2")
 
    f_ex = exact(zero(dt), mesh1,mesh2) 

    println("n=0 error=$(error1(f, f_ex))")

    for n=1:nt
#        println("dt=$dt")
        advection!(adv_x1, f,  v1, tan(dt/2))
        transpose!(ft, f)
        advection!(adv_x2, ft, v2, sin(dt))
        transpose!(f, ft)
        advection!(adv_x1, f, v1, tan(dt/2))
        f_ex = exact(n*dt, mesh1,mesh2) 
        println("n=$n error=$(error1(f, f_ex))")
        # if n%10 == 0
        #     savefile(abs.(f-f_ex),5,"save5_$(nt)_$n.png")
        #     savefile(abs.(f-f_ex),10,"save10_$(nt)_$n.png")
        #     savefile(abs.(f-f_ex),15,"save15_$(nt)_$n.png")
        #     savefile(abs.(f-f_ex),20,"save20_$(nt)_$n.png")
        # end
    end
 #   println(f)
    f

end
function rotation2d_2d(
    tf::T, 
    nt, 
    mesh1::UniformMesh{T}, 
    mesh2::UniformMesh{T}, 
    interp::InterpolationType
) where {T}

    dt = tf/nt

    n1 = mesh1.length
    x1min, x1max = mesh1.start, mesh1.stop
    delta1 = mesh1.step

    n2 = mesh2.length
    x2min, x2max = mesh2.start, mesh2.stop
    delta2 = mesh2.step

    fp  = zeros(T,(n1,n2))
    fp .= exact(zero(T), mesh1, mesh2)

#    println(f)

    
    adv = Advection2d( mesh1, mesh2, interp )

    tabv = [(-y,x) for y in mesh2.points, x in mesh1.points]
    
    for n=1:nt
        advection!(adv, fp, tabv, dt)
        println("n=$n error=$(error1(fp,exact(n*dt, mesh1,mesh2)))")
    end
 #   println(f)
    fp

end
# @testset "Rotation test with Lagrange2d advections " begin

#     tf, nt = 2π, 1000
    
#     mesh1 = UniformMesh(-pi, float(pi), 64; endpoint=false)
#     mesh2 = UniformMesh(-pi, float(pi), 64; endpoint=false)

#     @time lag= Lagrange2d(11)

#     println("norm lag = $(norm(lag.coef))")
    
#     @time fc = rotation2d_2d(tf, nt, mesh1, mesh2, lag)
#     fe = exact(tf, mesh1, mesh2)

 
#     err = error1(fc, fe)
#     println("err=$err")
#     @test err <  1e-1

# end

# @testset "Rotation test with Lagrange advections " begin
# function fct1()
#     tf, nt, nb = big(2π), 10, 256

#     println("nb=$nb")

#     mesh1 = UniformMesh(-big(5.0), big(5.0), nb; endpoint=false)
#     mesh2 = UniformMesh(-big(5.0), big(5.0), nb; endpoint=false)


#     @time lag= Lagrange(BigFloat, 51)

 
#     @time fc = rotation_2d(tf, nt, mesh1, mesh2, lag)
#     fe = exact(big"0.0", mesh1, mesh2)


#     err = error1(fc, fe)
#     println("err=$err")
#     @test err <  1e-1
# end

# fct1()

# end
# @testset "Rotation test with Lagrange advections " begin

# tf, nt = 2big(π), 100

# mesh1 = UniformMesh(-big(pi), big(pi), 128; endpoint=false)
# mesh2 = UniformMesh(-big(pi), big(pi), 128; endpoint=false)

# @time lag= Lagrange(21, iscirc=true)

# println("norm lag = $(norm(lag.coef))")
# println("len=$(size(mesh1.points,1)) nb=$nt order=$(size(lag.coef,1)-1)")

# @time fc = rotation_2d(tf, nt, mesh1, mesh2, lag)
# fe = exact(tf, mesh1, mesh2)


# err = error1(fc, fe)
# println("err=$err")
# @test err <  1e-1

# end
# @testset "Rotation test with Lagrange advections big" begin

#     tf, nt = 2big(π), 100
    
#     mesh1 = UniformMesh(-big(π), big(π), 64; endpoint=false)
#     mesh2 = UniformMesh(-big(π), big(π), 64; endpoint=false)

#     @time lag = Lagrange(BigFloat, 21)
    
#     println("norm lag = $(norm(lag.coef))")

#     @time fc = rotation_2d(tf, nt, mesh1, mesh2, lag )
#     fe = exact(tf, mesh1, mesh2)

#     err = error1(fc, fe)
#     println("err=$err")
#     @test err <  1e-1

# end

@testset "Rotation test with Lagrange advections " begin

tf, nt, nb = 2big(π), 10, 128

mesh1 = UniformMesh(-big(5.0), big(5.0), nb; endpoint=false)
mesh2 = UniformMesh(-big(5.0), big(5.0), nb; endpoint=false)

bsp = Lagrange(BigFloat, 31)

@time fc = rotation_2d(tf, nt, mesh1, mesh2, bsp)
fe = exact(tf, mesh1, mesh2)

err = error1(fc, fe)
println("err=$err")
@test err <  1e-7

end
@testset "Rotation test with BsplineLU advections" begin

    tf, nt, nb = 2big(π), 10, 128
    
    mesh1 = UniformMesh(-big(5.0), big(5.0), nb; endpoint=false)
    mesh2 = UniformMesh(-big(5.0), big(5.0), nb; endpoint=false)
    
    bsp = B_SplineLU(31, nb, BigFloat, iscirc=true)

    @time fc = rotation_2d(tf, nt, mesh1, mesh2, bsp)
    fe = exact(tf, mesh1, mesh2)

    err = error1(fc, fe)
    println("err=$err")
    @test err <  1e-10

end

@testset "Rotation test with BsplineFFT advections " begin

tf, nt, nb = 2big(π), 10, 128

mesh1 = UniformMesh(-big(5.0), big(5.0), nb; endpoint=false)
mesh2 = UniformMesh(-big(5.0), big(5.0), nb; endpoint=false)

bsp = B_SplineFFT(31, nb, BigFloat)

@time fc = rotation_2d(tf, nt, mesh1, mesh2, bsp)
fe = exact(tf, mesh1, mesh2)

err = error1(fc, fe)
println("err=$err")
@test err <  1e-10

end

# @testset "Rotation test with Bspline advections big" begin

#     tf, nt = 2big(π), 100
    
#     mesh1 = UniformMesh(-big(π), big(π), 128; endpoint=false)
#     mesh2 = UniformMesh(-big(π), big(π), 128; endpoint=false)
    
#     @time fc = rotation_2d(tf, nt, mesh1, mesh2, Bspline(5))
#     fe = exact(tf, mesh1, mesh2)

#     err = error1(fc, fe)
#     println("err=$err")
#     @test err <  1e-5

# end

