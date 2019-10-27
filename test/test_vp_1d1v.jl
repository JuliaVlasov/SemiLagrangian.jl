using SemiLagrangian

@testset " VP 1D1V " begin

eps    = 0.001
kx     = 0.5
degree = 5

xmin, xmax, nx =  0., 2π/kx, 32
vmin, vmax, nv = -6., 6., 64

mesh_x = UniformMesh( xmin, xmax, nx, endpoint=false )
mesh_v = UniformMesh( vmin, vmax, nv, endpoint=false )

advection_x! = PeriodicAdvection( mesh_x, Bspline(degree))
advection_v! = PeriodicAdvection( mesh_v, Bspline(degree))

ex  = zeros(Float64, nx)
rho = zeros(Float64, nx)

fxv = zeros(Float64, (nx, nv))
fvx = zeros(Float64, (nv, nx))

x = mesh_x.points
v = mesh_v.points

fxv = (1 .+ eps*cos.(kx .* x)) .* transpose(exp.(-.5*v.^2)) ./ sqrt(2π)

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
