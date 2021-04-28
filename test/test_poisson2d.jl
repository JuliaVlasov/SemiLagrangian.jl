using DoubleFloats
using LinearAlgebra

using SemiLagrangian:
    Advection,
    sizeall,
    AdvectionData,
    getdata,
    advection!,
    UniformMesh,
    AbstractInterpolation,
    Lagrange,
    B_SplineLU,
    B_SplineFFT,
    gettabmod,
    CachePrecal,
    interpolate!,
    compute_charge!,
    compute_elfield!,
    compute_ee,
    compute_ke
# """

#    exact(tf, mesh1, mesh2)

#    Julia function to compute exact solution "

# ```math
# \frac{d f}{dt} +  (y \frac{df}{delta1} - x \frac{df}{delta2}) = 0 
# ```

# """
# function exact!(f, mesh1::UniformMesh{T}, mesh2::UniformMesh{T}, tf::T) where {T}
#     for (i, x) in enumerate(mesh1.points), (j, y) in enumerate(mesh2.points)
#         s, c = sincos(tf)
#         xn, yn = c * x - s * y, + s * x + c * y
#         f[i,j] = exp(-5*((xn)^2+(yn+T(6//5))^2))
#     end
# end

# function test_rotation(
#     sz::NTuple{2,Int}, 
#     interp_sp::AbstractInterpolation{T}, 
#     interp_v::AbstractInterpolation{T},
#     nbdt::Int
# ) where {T}
#     spmin, spmax, nsp =  T(-5), T(5),  sz[1]
#     vmin, vmax, nv = -T(12//2), T(9//2), sz[2]

#     mesh_sp = UniformMesh( spmin, spmax, nsp)
#     mesh_v = UniformMesh( vmin, vmax, nv)

#     dt = T(2big(pi)/nbdt)
#     adv = Advection1d((mesh_sp,), (mesh_v,), (interp_sp,), (interp_v,), dt; tab_fct=[tan, sin, tan])
# #    adv = Advection1d((mesh_sp,), (mesh_v,), (interp_sp,), (interp_v,), dt; tab_coef=[1], tab_fct=[identity])
#     sz = sizeall(adv)
#     tabref = zeros(T,sz)
#     exact!(tabref, mesh_sp, mesh_v, T(0))

#     pvar = getrotationvar(adv)

#     advdata = Advection1dData(adv, tabref, pvar)

#     diffmax = 0
#     data = getdata(advdata)
#     for ind=1:nbdt
#         while advection!(advdata) end
#         exact!(tabref, mesh_sp, mesh_v, dt*ind)
#         diff = norm(data .- tabref, Inf)
#         diffmax = max(diffmax, diff)
#         @show ind, diff
#     end
#     println("test_rotation sz=$sz interp=$interp_sp, $interp_v nbdt=$nbdt diffmax=$diffmax")
#     return diffmax

# end

function test_poisson2d(
    sz::NTuple{2,Int},
    interp::Vector{I},
    t_max::T,
    nbdt::Int,
) where {T,I<:AbstractInterpolation{T}}
    spmin, spmax, nsp = T(0), 4T(pi), sz[1]
    vmin, vmax, nv = T(-6), T(6), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    dt = t_max / nbdt

    # dec1 = -2tan(dt/2) / step(mesh_sp) * mesh_v.points
    # dec2 = 2tan(dt/2) / step(mesh_v) * mesh_sp.points
    dec1 = zeros(T, sz)
    dec2 = zeros(T, sz)
    coef = 1 / sqrt(2T(pi))
    tabref = zeros(T, sz)
    for i = 1:sz[1], j = 1:sz[2]
        x = mesh_sp.points[i]
        y = mesh_v.points[j]
        dec1[i, j] = -dt / step(mesh_sp) * y
        tabref[i, j] = coef * exp(-0.5 * y^2) * (1 + 0.001 * cos(x / 2))
    end


    #    @show dec1, dec2

    cache = CachePrecal(interp, t_max)

    rho = zeros(T, sz[1])
    elf = zeros(T, sz[1])


    data = zeros(T, sz)

    copyto!(data, tabref)

    compute_charge!(rho, (mesh_sp,), data)
    compute_elfield!(elf, mesh_sp, rho)
    ee = compute_ee((mesh_sp,), (elf,))
    ke = compute_ke((mesh_sp,), (mesh_v,), data)

    println("0\t$ee\t$ke")


    buf = zeros(T, sz)
    diffmax = 0
    tabmod = gettabmod.(sz)
    for ind = 1:2nbdt
        dec2 = [dt / step(mesh_v) * elf[i] for i = 1:sz[1], j = 1:sz[2]]
        interpolate!(buf, data, ind->(dec1[ind], dec2[ind]), interp, tabmod=tabmod, cache=cache)
        copyto!(data, buf)

        compute_charge!(rho, (mesh_sp,), data)
        compute_elfield!(elf, mesh_sp, rho)

        ee = compute_ee((mesh_sp,), (elf,))
        ke = compute_ke((mesh_sp,), (mesh_v,), data)

        println("$(Float32(ind*dt))\t$ee\t$ke")
    end
    return diffmax
end

@testset "test swirling" begin
    T = Double64
    @time @test test_poisson2d((256, 200), [Lagrange(11, T),Lagrange(11, T)] , T(10), 10) < 1e-3
end
