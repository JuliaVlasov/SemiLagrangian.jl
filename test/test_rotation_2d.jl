abstract type InterpolationType end
abstract type AbstractAdvection end
include("../src/mesh.jl")
include("../src/bspline_periodic.jl")
include("../src/advection.jl")
include("../src/lagrange.jl")





using LinearAlgebra

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
        xn = cos(tf) * x + sin(tf) * y
        yn = - sin(tf) * x + cos(tf) * y
        f[i,j] = (xn-0.3)^2+(yn+0.5)^2 < 0.03 ? 1.0 : 0.0
#        f[i,j] = exp(-(xn-0.3)*(xn-0.3)/0.2)*exp(-(yn+0.5)*(yn+0.5)/0.2)
#       f[i,j] = exp(-(xn-1)*(xn-1)*10)*exp(-(yn)*(yn)*10)
#       f[i,j] = exp(-(sin(xn)+0.4)^2)*exp(-(cos(yn)-0.5)^2)
#       f[i,j] = exp(-sin(xn-0.3)*sin(xn-0.3)/0.6)*exp(-cos(yn+0.5)*cos(yn+0.5)/0.6)
    end

    f

end

" Function to compute error "
function error1(f, f_exact)
    maximum(abs.(f .- f_exact))
end

function rotation_2d(
    tf, 
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

    f  = zeros(T,(n1,n2))
    f .= exact(zero(T), mesh1, mesh2)
    ft = zeros(T,(n2,n1))
    transpose!(ft, f)

#    println(f)

    
    adv_x1 = Advection( mesh1, interp )
    adv_x2 = Advection( mesh2, interp )

    v1 = - collect(mesh2.points)
    v2 = + collect(mesh1.points)
    
    for n=1:nt
        advection!(adv_x1, f,  v1, tan(dt/2))
        transpose!(ft, f)
        advection!(adv_x2, ft, v2, sin(dt))
        transpose!(f, ft)
        advection!(adv_x1, f, v1, tan(dt/2))
        println("n=$n error=$(error1(f,exact(n*dt, mesh1,mesh2)))")
    end
 #   println(f)
    f

end

@testset "Rotation test with LagrangeNew advections " begin

    tf, nt = 2π, 100
    
    mesh1 = UniformMesh(-pi, float(pi), 64; endpoint=false)
    mesh2 = UniformMesh(-pi, float(pi), 64; endpoint=false)

    @time lag= LagrangeNew(30, granularity=8)

    println("norm lag = $(norm(lag.coef))")
    
    @time fc = rotation_2d(tf, nt, mesh1, mesh2, lag)
    fe = exact(tf, mesh1, mesh2)

 
    err = error1(fc, fe)
    println("err=$err")
    @test err <  1e-1

end
# @testset "Rotation test with LagrangeNew advections big" begin

#     tf, nt = 2big(π), 100
    
#     mesh1 = UniformMesh(-big(π), big(π), 64; endpoint=false)
#     mesh2 = UniformMesh(-big(π), big(π), 64; endpoint=false)

#     @time lag = LagrangeNew(BigFloat, 40, granularity=10)
    
#     println("norm lag = $(norm(lag.coef))")

#     @time fc = rotation_2d(tf, nt, mesh1, mesh2, lag )
#     fe = exact(tf, mesh1, mesh2)

#     err = error1(fc, fe)
#     println("err=$err")
#     @test err <  1e-1

# end
# @testset "Rotation test with Bspline advections " begin

#     tf, nt = 2π, 100
    
#     mesh1 = UniformMesh(-π, float(π), 128; endpoint=false)
#     mesh2 = UniformMesh(-π, float(π), 128; endpoint=false)
    
#     @time fc = rotation_2d(tf, nt, mesh1, mesh2, Bspline(5))
#     fe = exact(tf, mesh1, mesh2)

#     err = error1(fc, fe)
#     println("err=$err")
#     @test err <  1e-4

# end

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

