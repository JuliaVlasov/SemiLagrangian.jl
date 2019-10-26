# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.2.0
#     language: julia
#     name: julia-1.2
# ---

using Plots, ProgressMeter

# +
abstract type InterpolationType end
abstract type AbstractAdvection end

include("../src/mesh.jl")
include("../src/geometry.jl")
include("../src/bspline_periodic.jl")
include("../src/splinenn.jl")
include("../src/splinepp.jl")
include("../src/banded_matrix.jl")
include("../src/spline_1d.jl")
include("../src/spline_interpolator_1d.jl")
include("../src/advection.jl")

# +
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
        f[i,j] = exp(-(xn-1)*(xn-1)/0.1)*exp(-(yn-1)*(yn-1)/0.1)
    end

    f

end
# -

" Function to compute error "
function error1(f, f_exact)
    maximum(abs.(f .- f_exact))
end

# +
tf, nt = 2π, 100

mesh1 = UniformMesh(-π, π, 256; endpoint=false)
mesh2 = UniformMesh(-π, π, 128; endpoint=false)

n1 = mesh1.length
n2 = mesh2.length

f  = zeros(Float64,(n1,n2))
f .= exact(0.0, mesh1, mesh2)
ft = zeros(Float64,(n2,n1))
transpose!(ft, f)
advection_x1! = PeriodicAdvection( mesh1, Bspline(5) )
advection_x2! = PeriodicAdvection( mesh2, Bspline(5) )

v1 = collect(mesh2.points)
v2 = - collect(mesh1.points)
dt = tf / nt

 @gif for it = 1:nt
    
    advection_x1!( f,  v1, tan(dt/2))
    transpose!(ft, f)
    advection_x2!( ft,  v2, sin(dt))
    transpose!(f, ft)
    advection_x1!( f,  v1, tan(dt/2))
    p = contour(mesh2.points, mesh1.points, f)
    plot!(p[1]; clims=(0.,1.), aspect_ratio=:equal)
    plot!( sqrt(2) .* cos.(-pi:0.1:pi),  sqrt(2) .* sin.(-pi:0.1:pi))
    
end
