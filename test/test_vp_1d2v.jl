
import LinearAlgebra: transpose!

@testset " VP 1D2V " begin

function transpose!(fvx::Array{Float64,3}, fxv::Array{Float64,3})

    permutedims!(fvx, fxv, [2, 3, 1])

end

function transpose!(fxv::Array{Float64,3}, fvx::Array{Float64,3})

    permutedims!(fxv, fvx, [3, 1, 2])

end


eps    = 0.001
kx     = 0.5
degree = 5

x1min, x1max, nx1 =  0., 2π/kx, 32
v1min, v1max, nv1 = -6., 6., 64
v2min, v2max, nv2 = -6., 6., 64

mesh_x = UniformMesh( xmin, xmax, nx, endpoint=false )
mesh_v = RectMesh2D( v1min, v1max, nv1, v2min, v2max, nv2)

advection_x! = PeriodicAdvection( mesh_x, BsplineOld(degree))
advection_v! = Advection( mesh_v, BsplineOld(degree))

ex  = zeros(Float64, nx)
rho = zeros(Float64, nx)

fxv = zeros(Float64, (nx1, nv1, nv2))
fvx = zeros(Float64, (nv1, nv2, nx1))

x1 = mesh_x.points
v1, v2 = mesh_v.points

for i in eachindex(x1), j in eachindex(v1), k in eachindex(v2)
   fxv[i,j,k] = (1 + eps*cos(kx * x1[i])) * exp(-.5*(v1[j]^2+v2[k]^2)) / sqrt(2π)
end

transpose!(fvx, fxv)

tspan  = LinRange(0, 60, 600)
dt = 0.1

compute_charge!(rho, mesh_v, fvx)

compute_e!(ex, mesh_x, rho)

for t in tspan

  advection_x!(fxv, v, 0.5dt)

  transpose!(fvx, fxv)

  compute_charge!(rho, mesh_v, fvx)
     
  compute_e!(ex, mesh_x, rho)

  advection_v!(fvx, ex, dt)

  transpose!(fxv, fvx)

  advection_x!(fxv, v, 0.5dt)

end 

@test true

end
