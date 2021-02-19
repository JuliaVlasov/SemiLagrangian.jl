# Quickstart

## Rotation
First simulation using semi-lagrangian method to get a rotation

```julia

using SemiLagrangian
using LinearAlgebra

function exact!(f, mesh1::UniformMesh{T}, mesh2::UniformMesh{T}, tf::T) where {T}
    for (i, x) in enumerate(mesh1.points), (j, y) in enumerate(mesh2.points)
        s, c = sincos(tf)
        xn, yn = c * x - s * y, s * x + c * y
        f[i,j] = exp(-13*((xn)^2+(yn+T(6//5))^2))
    end
end

mesh_sp = UniformMesh( -5.0, 5.0, 256)
mesh_v = UniformMesh( -5.0, 5.0, 256)
nbdt=50
dt = 2pi/nbdt
interp_sp = Lagrange(11)
interp_v = Lagrange(11)

adv = Advection((mesh_sp,), (mesh_v,), (interp_sp,), (interp_v,), dt; tab_fct=[tan, sin, tan])
sz = sizeall(adv)
tabref = zeros(sz)
exact!(tabref, mesh_sp, mesh_v, 0.0)

pvar = getrotationvar(adv)

advdata = AdvectionData(adv, tabref, pvar)

diffmax = 0
data = getdata(advdata)
for ind=1:nbdt
    while advection!(advdata) end
    exact!(tabref, mesh_sp, mesh_v, dt*ind)
    diff = norm(data .- tabref, Inf)
    diffmax = max(diffmax, diff)
    println("ind=$ind sz=$sz interp=$interp_sp, $interp_v nbdt=$nbdt diff,diffmax=$diff,$diffmax")
end

```

## Vlasov/poisson

For Vlasov Poisson equation 
```julia
using SemiLagrangian

mesh_sp = UniformMesh( 0.0, 4pi, 32)
mesh_v = UniformMesh( -6.0, 6.0, 32)

nbdt=50
dt = 0.1
epsilon = 0.5

interp = Lagrange(7)

adv = Advection( (mesh_sp, mesh_sp), (mesh_v, mesh_v), (interp, interp), (interp, interp), dt)

fct_sp(x)=epsilon * cos(x/2) + 1
fct_v(v)=exp( - v^2 / 2)/sqrt(2pi)
lgn_sp = fct_sp.(mesh_sp.points)
lgn_v = fct_v.(mesh_v.points)

data = dotprod((lgn_sp, lgn_sp, lgn_v, lgn_v))

pvar = getpoissonvar(adv)

advdata = AdvectionData(adv, data, pvar)

t=0
compute_charge!(advdata)
compute_elfield!(advdata)
ee = compute_ee(advdata)
ke = compute_ke(advdata)
println("$t\t$ee\t$ke")

for ind=1:nbdt
    while advection!(advdata) end
    t = Float32(dt*ind)
    compute_charge!(advdata)
    compute_elfield!(advdata)
    ee = compute_ee(advdata)
    ke = compute_ke(advdata)
    println("$t\t$ee\t$ke")
end

```

