using DoubleFloats
using LinearAlgebra

using SemiLagrangian: Advection, sizeall, AdvectionData, getdata, advection!, UniformMesh, getrotationvar, AbstractInterpolation, 
Lagrange, B_SplineLU, B_SplineFFT
"""

   exact(tf, mesh1, mesh2)

   Julia function to compute exact solution "

```math
\frac{d f}{dt} +  (y \frac{df}{delta1} - x \frac{df}{delta2}) = 0 
```

"""
function exact!(f, mesh1::UniformMesh{T}, mesh2::UniformMesh{T}, tf::T) where {T}
    for (i, x) in enumerate(mesh1.points), (j, y) in enumerate(mesh2.points)
        s, c = sincos(tf)
        xn, yn = c * x - s * y, s * x + c * y
        f[i,j] = exp(-13*((xn)^2+(yn+T(6//5))^2))
    end
end

function test_rotation(
    sz::NTuple{2,Int}, 
    interp_sp::AbstractInterpolation{T}, 
    interp_v::AbstractInterpolation{T},
    nbdt::Int
) where {T}
    spmin, spmax, nsp =  T(-5), T(5),  sz[1]
    vmin, vmax, nv = -T(5), T(5), sz[2]

    mesh_sp = UniformMesh( spmin, spmax, nsp)
    mesh_v = UniformMesh( vmin, vmax, nv)

    dt = T(2big(pi)/nbdt)
    adv = Advection((mesh_sp,), (mesh_v,), (interp_sp,), (interp_v,), dt; tab_fct=[tan, sin, tan])
    sz = sizeall(adv)
    tabref = zeros(T,sz)
    exact!(tabref, mesh_sp, mesh_v, T(0))

    pvar = getrotationvar(adv)

    advdata = AdvectionData(adv, tabref, pvar)

    diffmax = 0
    data = getdata(advdata)
    for ind=1:nbdt
        while advection!(advdata) end
        exact!(tabref, mesh_sp, mesh_v, dt*ind)
        diff = norm(data .- tabref, Inf)
        diffmax = max(diffmax, diff)
    end
    println("test_rotation sz=$sz interp=$interp_sp, $interp_v nbdt=$nbdt diffmax=$diffmax")
    return diffmax

end


@testset "test rotation" begin
    T = Float64
    @time @test test_rotation((200, 200), Lagrange(5, T), Lagrange(5, T), 11) < 1e-3
    @time @test test_rotation((128, 128), B_SplineLU(5, 128, T), B_SplineLU(5, 128, T), 11) < 1e-3
    @time @test test_rotation((128, 128), B_SplineFFT(5, 128, T), B_SplineFFT(5, 128, T), 11) < 1e-3
    T = Double64
    @time @test test_rotation((200, 200), Lagrange(15, T), Lagrange(15, T), 11) < 1e-7
    @time @test test_rotation((128, 128), B_SplineLU(15, 128, T), B_SplineLU(15, 128, T), 11) < 1e-8
    @time @test test_rotation((128, 128), B_SplineFFT(15, 128, T), B_SplineFFT(15, 128, T), 11) < 1e-8
    T = BigFloat
    @time @test test_rotation((200, 200), Lagrange(25, T), Lagrange(25, T), 11) < 1e-10
    @time @test test_rotation((128, 128), B_SplineLU(25, 128, T), B_SplineLU(25, 128, T), 11) < 1e-11
    @time @test test_rotation((128, 128), B_SplineFFT(25, 128, T), B_SplineFFT(25, 128, T), 11) < 1e-11

end
