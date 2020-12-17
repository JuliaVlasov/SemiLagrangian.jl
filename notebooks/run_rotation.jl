# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Julia 1.4.1
#     language: julia
#     name: julia-1.4
# ---

using LinearAlgebra
# using ProgressMeter
# using Plots

# +
import Base.Threads: @spawn, @sync, nthreads, threadid
# include("../src/mpiinterface.jl")
include("../src/nompiinterface.jl")
include("../src/advection.jl")
include("../src/rotation.jl")
include("../src/spline.jl")
include("../src/bspline.jl")
include("../src/bsplinelu.jl")
include("../src/bsplinefft.jl")
include("../src/lagrange.jl")
include("../src/interpolation.jl")

using DoubleFloats

function printout(advd::AdvectionData{T,Nsp,Nv,Nsum,timeopt}, str) where {T,Nsp,Nv,Nsum,timeopt}
    if timeopt != MPIOpt || advd.adv.mpid.ind == 1
        println(str)
    end
end

"""

   exact(tf, mesh1, mesh2)

   Julia function to compute exact solution "

```math
\frac{d f}{dt} +  (y \frac{df}{delta1} - x \frac{df}{delta2}) = 0 
```

"""
function exact!(f, mesh1::UniformMesh{T}, mesh2::UniformMesh{T}, tf::T) where {T}

    for (i, x) in enumerate(mesh1.points), (j, y) in enumerate(mesh2.points)
        xn = cos(tf) * x - sin(tf) * y
        yn = sin(tf) * x + cos(tf) * y
        f[i,j] = exp(-13*((xn)^2+(yn+big"1.2")^2))
    end
    f
end

function trace_diffrotation(advd::AdvectionData{T,Nsp,Nv,Nsum,timeopt}, t::T) where{T,Nsp,Nv,Nsum,timeopt}

    if t == 0
        printout(advd, "#time\tnorm_2\tnorm_inf")
    end
    ptr = pointer(advd.bufdata)
    f = unsafe_wrap(Array, ptr, sizeall(advd.adv))
    exact!(f, advd.adv.t_mesh_sp[1], advd.adv.t_mesh_v[1], t)
    n2 = Float64(norm(advd.data-f))
    ninf = Float64(norm(advd.data-f, Inf))
    printout(advd, "$(Float32(t))\t$n2\t$ninf")

end

function rotation(advd::AdvectionData, nbdt)

    # global cl_obs
    # clockreset(cl_obs)

    dt = advd.adv.dt_base
    trace_diffrotation(advd, zero(dt))
    for i=1:nbdt
        while advection!(advd)
        end
        trace_diffrotation(advd, i*dt)
    end
    println("#  end")
# printall(cl_obs)
end
function rotation1_1(T::DataType, nbdt, timeopt; sz=(64,64), interp=Lagrange(T, 21))
    epsilon = T(0.5)
    dt = T(2big(pi)/nbdt)

    spmin, spmax, nsp =  T(-5), T(5),  sz[1]
    vmin, vmax, nv = -T(5), T(5), sz[2]

    mesh_sp = UniformMesh( spmin, spmax, nsp, endpoint = false)
    mesh_v = UniformMesh( vmin, vmax, nv, endpoint = false )
    
    adv = Advection((mesh_sp,), (mesh_v,), (interp,), (interp,), dt, timeopt=timeopt, tab_fct=[tan,sin,tan])

    data = zeros(T, sizeall(adv))

    exact!(data, mesh_sp, mesh_v, zero(T))

    pvar = getrotationvar(adv)

    advd = AdvectionData(adv, data, pvar)

    printout(advd, "# dt=$(Float64(dt)) eps=$(Float64(epsilon)) size_x=$nsp size_v=$nv")
    printout(advd, "# sp : from $(Float64(mesh_sp.start)) to $(Float64(mesh_sp.stop))")
    printout(advd, "# v : from $(Float64(mesh_v.start)) to $(Float64(mesh_v.stop))")
    printout(advd, "# interpolation : $(get_type(interp)) order=$(get_order(interp))")
    printout(advd, "# type=$T precision = $(precision(T))")
    printout(advd, "# timeopt=$timeopt")
    if timeopt == SimpleThreadsOpt || timeopt == SplitThreadsOpt
        printout(advd, "# nb threads : $(Threads.nthreads())")
    elseif timeopt == MPIOpt
        printout(advd,"# nb process : $(adv.mpid.nb)")
    else
        printout(advd, "# monothread version")
    end
    printout(advd, "typeof(data)=$(typeof(data)) size(data)=$(size(data))")

    # advdata = AdvectionData(adv, data, pvar)

    rotation(advd, nbdt)
end   
# function landau2_2(T::DataType, nbdt, timeopt; sz=(32,32,32,32), dt = big"0.1", interp=Lagrange(T, 31))
#     epsilon = T(0.5)
#     dt = T(dt)

#     sp1min, sp1max, nsp1 =  T(0), T(4big(pi)),  sz[1]
#     v1min, v1max, nv1 = -T(6.), T(6.), sz[3]

#     mesh1_sp = UniformMesh( sp1min, sp1max, nsp1, endpoint = false)
#     mesh1_v = UniformMesh( v1min, v1max, nv1, endpoint = false )

#     sp2min, sp2max, nsp2 =  T(0), T(4big(pi)),  sz[2]
#     v2min, v2max, nv2 = -T(6.), T(6.), sz[4]

#     mesh2_sp = UniformMesh( sp2min, sp2max, nsp2, endpoint = false)
#     mesh2_v = UniformMesh( v2min, v2max, nv2, endpoint = false )

 
#     adv = Advection((mesh1_sp, mesh2_sp), (mesh1_v, mesh2_v), (interp,interp,), (interp,interp,), dt, timeopt=timeopt)

#     fct_sp(x)=epsilon * cos(x/2) + 1
#     fct_v(v)=exp( - v^2 / 2)/sqrt(2T(pi))
#     lgn1_sp = fct_sp.(mesh1_sp.points)
#     lgn1_v = fct_v.(mesh1_v.points)
#     lgn2_sp = fct_sp.(mesh2_sp.points)
#     lgn2_v = fct_v.(mesh2_v.points)

#     data = dotprod((lgn1_sp, lgn2_sp, lgn1_v, lgn2_v))
    

#     pvar = getpoissonvar(adv)

#     advd = AdvectionData(adv, data, pvar)
#     # advdata = AdvectionData(adv, data, pvar)
#     printout(advd, "# dt=$(Float32(dt)) eps=$(Float64(epsilon)) size1_sp=$nsp1 size2_sp=$nsp2 size_v1=$nv1 size_v2=$nv2")
#     printout(advd, "# sp1 : from $(Float64(mesh1_sp.start)) to $(Float64(mesh1_sp.stop))")
#     printout(advd, "# sp2 : from $(Float64(mesh2_sp.start)) to $(Float64(mesh2_sp.stop))")
#     printout(advd, "# v1 : from $(Float64(mesh1_v.start)) to $(Float64(mesh1_v.stop))")
#     printout(advd, "# v2 : from $(Float64(mesh2_v.start)) to $(Float64(mesh2_v.stop))")
#     printout(advd, "# interpolation : $(get_type(interp)) order=$(get_order(interp))")
#     printout(advd, "# type=$T precision = $(precision(T))")
#     printout(advd, "# timeopt=$timeopt")
#     if timeopt == SimpleThreadsOpt || timeopt == SplitThreadsOpt
#         printout(advd, "# nb threads : $(Threads.nthreads())")
#     elseif timeopt == MPIOpt
#         printout(advd,"# nb process : $(adv.mpid.nb)")
#     else
#         printout(advd, "# monothread version")
#     end
#     printout(advd, "typeof(data)=$(typeof(data)) size(data)=$(size(data))")

#     landau(advd, nbdt)
# end

T=Double64
nbdt=1000
dt=T(2big(pi))/nbdt
tab_fct=[tan,sin,tan]
# landau1_1(T, 50, NoTimeOpt, sz=(64,128))
rotation1_1(T, nbdt, NoTimeOpt, sz=(256,256), interp=Lagrange(T, 101))
