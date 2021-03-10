using DoubleFloats
using LinearAlgebra

using SemiLagrangian: Advection, sizeall, AdvectionData, getdata, advection!, UniformMesh, gettranslationvar, AbstractInterpolation, interpolate!,
Lagrange, B_SplineLU, B_SplineFFT, Lagrange2d
"""

   exact(tf, mesh1, mesh2)

   Julia function to compute exact solution "

```math
\frac{d f}{dt} +  (y \frac{df}{delta1} - x \frac{df}{delta2}) = 0 
```

"""
function exact!(f, decall::NTuple{Nsum,T}, t::T) where {T, Nsum}
    for ind in CartesianIndices(f)
       f[ind] = exp( -sum( [sin( 2T(pi)*(ind.I[i]+decall[i]*t)/size(f,i) +T(big"0.56")*i^2) for i=1:Nsum] ) )
    end
end

function test_translation(
    sz::NTuple{2,Int},
    interp_sp::AbstractInterpolation{T}, 
    interp_v::AbstractInterpolation{T},
    nbdt::Int
) where {T}
    spmin, spmax, nsp =  T(-5), T(5),  sz[1]
    vmin, vmax, nv = -T(5.55), T(5), sz[2]

    mesh_sp = UniformMesh( spmin, spmax, nsp)
    mesh_v = UniformMesh( vmin, vmax, nv)

    dt = T(1)
    adv = Advection((mesh_sp,), (mesh_v,), (interp_sp,), (interp_v,), dt)
    tabref = zeros(T,sz)
 
    v1=T(big"0.83545655467782872872782870029282982828737872878776717190927267611111")
    v2=T(-big"0.945678101929276765616767176761671771766717828781828998101092981877817176")

    exact!(tabref, (v1, v2), T(0))

    pvar = gettranslationvar((v1, v2))

    advdata = AdvectionData(adv, tabref, pvar)

    diffmax = 0
    data = getdata(advdata)
    for ind=1:nbdt
        while advection!(advdata) end
        exact!(tabref, (v1, v2), T(ind))
        diff = norm(data .- tabref, Inf)
#        @show v1, v2, ind, diff
        diffmax = max(diffmax, diff)
    end
    println("test_translation sz=$sz interp=$interp_sp, $interp_v nbdt=$nbdt diffmax=$diffmax")
    return diffmax

end
function test_translation(
    sz::NTuple{2,Int},
    interp::Lagrange2d{T}, 
    nbdt::Int
) where {T}
    spmin, spmax, nsp =  T(-5), T(5),  sz[1]
    vmin, vmax, nv = -T(6.5), T(5), sz[2]

    mesh_sp = UniformMesh( spmin, spmax, nsp)
    mesh_v = UniformMesh( vmin, vmax, nv)

    dt = T(1)
 
    v1=T(big"5.82545655467782872872782870029282982828737872878776717190927267611111")
    v2=T(-big"3.915678101929276765616767176761671771766717828781828998101092981877817176")

    @show v1, v2

    tabref = zeros(T,sz)
    exact!(tabref, (v1, v2), T(0))

    data = copy(tabref)

    buf = zeros(T,sz)

    dec1=fill(v1, sz)
    dec2=fill(v2, sz)
    diffmax=0
    for ind=1:nbdt
        interpolate!(buf, data, (dec1, dec2), interp)
        copyto!(data, buf)
        exact!(tabref, (v1, v2), T(ind))
        diff = norm(data .- tabref, Inf)
        diffmax = max(diffmax, diff)
        @show ind, diff
    end
    println("test_translation sz=$sz interp=$interp nbdt=$nbdt diffmax=$diffmax")
    return diffmax

end


@testset "test translation" begin
    T = Float64
    @time @test test_translation((200, 220), Lagrange2d(5, T), 11) < 1e-7
    @time @test test_translation((200, 220), Lagrange(5, T), Lagrange(5, T), 11) < 1e-7
    @time @test test_translation((128, 256), B_SplineLU(5, 128, T), B_SplineLU(5, 256, T), 11) < 1e-8
    @time @test test_translation((128, 256), B_SplineFFT(5, 128, T), B_SplineFFT(5, 256, T), 11) < 1e-8
    T = Double64
    @time @test test_translation((200, 220), Lagrange2d(15, T), 11) < 1e-18
    @time @test test_translation((200, 220), Lagrange(15, T), Lagrange(15, T), 11) < 1e-18
    @time @test test_translation((128, 256), B_SplineLU(15, 128, T), B_SplineLU(15, 256, T), 11) < 1e-22
    @time @test test_translation((128, 256), B_SplineFFT(15, 128, T), B_SplineFFT(15, 256, T), 11) < 1e-22
    T = BigFloat
    @time @test test_translation((200, 190), Lagrange2d(21, T), 11) < 1e-23
    @time @test test_translation((200, 190), Lagrange(25, T), Lagrange(25, T), 11) < 1e-28
    @time @test test_translation((128, 256), B_SplineLU(25, 128, T), B_SplineLU(25, 256, T), 11) < 1e-34
    @time @test test_translation((128, 256), B_SplineFFT(25, 128, T), B_SplineFFT(25, 256, T), 11) < 1e-34
end
