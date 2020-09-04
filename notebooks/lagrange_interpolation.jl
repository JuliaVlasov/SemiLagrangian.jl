# -*- coding: utf-8 -*-
using Plots
using ProgressMeter
include("../src/lagrange_interpolation.jl")

nx, nv = 101, 101
xmin

f( x ) =  exp( - x^2 / 2)


# +
xp = LinRange(a, b, n+1)[1:end-1] # remove end point
dx = xp[2] - xp[1]
yp = f.(xp) ;
xs = collect(xp)
ys = similar(yp)

pm = Progress(1000)
@gif for i in 1:1000
    alpha = -0.1
    interpolate!(ys, yp, alpha, LagrangeOld(5));
    yp .= ys
    xs .= a .+ mod.(xs .- a .+ alpha * dx, b-a)
    plot(xp, f.(xs), label="true") 
    plot!(xp, yp, label="computed")
    ylims!(0,1)
    next!(pm)
end
    

# -


