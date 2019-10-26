using SemiLagrangian


"""

   exact(tf, mesh1, mesh2)

   Julia function to compute exact solution "

```math
\frac{d f}{dt} +  (y \frac{df}{delta1} - x \frac{df}{delta2}) = 0 
```

"""
function exact(tf, mesh1::UniformMesh, mesh2::UniformMesh)

    f = zeros(Float64,(mesh1.length,mesh2.length))
    for (i, x) in enumerate(mesh1.points), (j, y) in enumerate(mesh2.points)
        xn = cos(tf) * x - sin(tf) * y
        yn = sin(tf) * x + cos(tf) * y
        f[i,j] = exp(-(xn-1)*(xn-1)/0.1)*exp(-(yn)*(yn)/0.1)
    end

    f

end

" Function to compute error "
function error1(f, f_exact)
    maximum(abs.(f .- f_exact))
end

function rotation_2d(tf, nt, mesh1::UniformMesh, mesh2::UniformMesh)

    dt = tf/nt

    n1 = mesh1.length
    x1min, x1max = mesh1.start, mesh1.stop
    delta1 = mesh1.step

    n2 = mesh2.length
    x2min, x2max = mesh2.start, mesh2.stop
    delta2 = mesh2.step

    f  = zeros(Float64,(n1,n2))
    f .= exact(0.0, mesh1, mesh2)
    ft = zeros(Float64,(n2,n1))
    transpose!(ft, f)
    
    advection_x1! = PeriodicAdvection( mesh1, Bspline(5) )
    advection_x2! = PeriodicAdvection( mesh2, Bspline(5) )

    v1 = - collect(mesh2.points)
    v2 = + collect(mesh1.points)
    
    for n=1:nt

        advection_x1!( f,  v1, tan(dt/2))
        transpose!(ft, f)
        advection_x2!( ft, v2, sin(dt))
        transpose!(f, ft)
        advection_x1!( f, v1, tan(dt/2))

    end

    f

end

@testset "Rotation test with Bspline advections " begin

    tf, nt = 2π, 100
    
    mesh1 = UniformMesh(-π, π, 256; endpoint=false)
    mesh2 = UniformMesh(-π, π, 256; endpoint=false)
    
    @time fc = rotation_2d(tf, nt, mesh1, mesh2)
    fe = exact(tf, mesh1, mesh2)

    @test error1(fc, fe) <  1e-5

end
