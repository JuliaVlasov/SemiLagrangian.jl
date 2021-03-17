using DoubleFloats
using LinearAlgebra

using SemiLagrangian: Advection, sizeall, AdvectionData, getdata, advection!, UniformMesh, getrotationvar, AbstractInterpolation, 
Lagrange, B_SplineLU, B_SplineFFT, interpolate!
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
        xn, yn = c * x - s * y, + s * x + c * y
        f[i,j] = exp(-5*((xn)^2+(yn+T(6//5))^2))
    end
end

function test_rotation(
    sz::NTuple{2,Int}, 
    interp_sp::AbstractInterpolation{T}, 
    interp_v::AbstractInterpolation{T},
    nbdt::Int
) where {T}
    spmin, spmax, nsp =  T(-5), T(5),  sz[1]
    vmin, vmax, nv = -T(12//2), T(9//2), sz[2]

    mesh_sp = UniformMesh( spmin, spmax, nsp)
    mesh_v = UniformMesh( vmin, vmax, nv)

    dt = T(2big(pi)/nbdt)
    adv = Advection((mesh_sp,), (mesh_v,), (interp_sp,), (interp_v,), dt; tab_fct=[tan, sin, tan])
#    adv = Advection((mesh_sp,), (mesh_v,), (interp_sp,), (interp_v,), dt; tab_coef=[1], tab_fct=[identity])
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
        @show ind, diff
    end
    println("test_rotation sz=$sz interp=$interp_sp, $interp_v nbdt=$nbdt diffmax=$diffmax")
    return diffmax

end

function test_rotation(
    sz::NTuple{2,Int}, 
    interp::AbstractInterpolation{T, edge, order, 2}, 
    nbdt::Int
) where{T, edge, order}
    spmin, spmax, nsp =  T(-4.5), T(5),  sz[1]
    vmin, vmax, nv = -T(5.2), T(5), sz[2]

    mesh_sp = UniformMesh( spmin, spmax, nsp)
    mesh_v = UniformMesh( vmin, vmax, nv)

    dt = T(2big(pi)/nbdt)

    # dec1 = -2tan(dt/2) / step(mesh_sp) * mesh_v.points
    # dec2 = 2tan(dt/2) / step(mesh_v) * mesh_sp.points
    dec1 = zeros(T, sz)
    dec2 = zeros(T, sz)
    tgdt = 2tan(dt/2)
    coef = 1/(1+tgdt^2/4)
    for i=1:sz[1], j=1:sz[2]
        x = mesh_sp.points[i] 
        y = mesh_v.points[j]
        v = [- tgdt * coef * (y + tgdt*x/2)/ step(mesh_sp) , (tgdt * coef * (x - tgdt*y/2))/step(mesh_v)]
        # nv = norm(v)
        # nr = norm([x,y])
        # if nv != 0
        #     v *= cor*nr/nv
        # end
        dec1[i,j], dec2[i,j] = v
    end

    @show minimum(dec1), maximum(dec1)
    @show minimum(dec2), maximum(dec2)

#    @show dec1, dec2

    tabref = zeros(T, sz)
    data = zeros(T, sz)

    exact!(tabref, mesh_sp, mesh_v, T(0))
    copyto!(data, tabref)

    buf = zeros(T, sz)
    diffmax=0
    for ind=1:nbdt
        interpolate!(buf, data, (x -> dec1[x],y -> dec2[y]), interp)
        copyto!(data, buf)
        exact!(tabref, mesh_sp, mesh_v, dt*ind)
        diff = norm(data .- tabref, Inf)
        diffmax = max(diffmax, diff)
        @show ind, diff

    end
    println("test_rotation sz=$sz interp=$interp nbdt=$nbdt diffmax=$diffmax")
    return diffmax
end

@testset "test rotation" begin
    T = Float64
    @time @test test_rotation((100, 122), Lagrange(5, T, nd=2), 11) < 1e-3
    @time @test test_rotation((200, 222), Lagrange(5, T), Lagrange(5, T), 11) < 1e-3
    @time @test test_rotation((128, 256), B_SplineLU(5, 128, T), B_SplineLU(5, 256, T), 11) < 1e-3
    @time @test test_rotation((128, 256), B_SplineFFT(5, 128, T), B_SplineFFT(5, 256, T), 11) < 1e-3
    T = Double64
    @time @test test_rotation((100, 122), Lagrange(15, T, nd=2), 11) < 1e-6
    @time @test test_rotation((200, 220), Lagrange(15, T), Lagrange(15, T), 11) < 1e-7
    @time @test test_rotation((128, 256), B_SplineLU(15, 128, T), B_SplineLU(15, 256, T), 11) < 1e-8
    @time @test test_rotation((128, 256), B_SplineFFT(15, 128, T), B_SplineFFT(15, 256, T), 11) < 1e-8
    T = BigFloat
    @time @test test_rotation((100, 122), Lagrange(21, T, nd=2), 11) < 1e-7
    @time @test test_rotation((200, 220), Lagrange(25, T), Lagrange(25, T), 11) < 1e-10
    @time @test test_rotation((128, 256), B_SplineLU(25, 128, T), B_SplineLU(25, 256, T), 11) < 1e-11
    @time @test test_rotation((128, 256), B_SplineFFT(25, 128, T), B_SplineFFT(25, 256, T), 11) < 1e-11
end
