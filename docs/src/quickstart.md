# Quickstart

## Rotation
First simulation using semi-lagrangian method to get a rotation

```@example rotation

import SemiLagrangian: magicsplit

using SemiLagrangian
using LinearAlgebra

function exact!(f, mesh1::UniformMesh{T}, mesh2::UniformMesh{T}, tf::T) where {T}
    for (i, x) in enumerate(mesh1.points), (j, y) in enumerate(mesh2.points)
        s, c = sincos(tf)
        xn, yn = c * x - s * y, s * x + c * y
        f[i,j] = exp(-13*((xn)^2+(yn+T(6//5))^2))
    end
end

function run_rotation()

    mesh_sp = UniformMesh( -5.0, 5.0, 256)
    mesh_v = UniformMesh( -5.0, 5.0, 256)
    nbdt=50
    dt = 2pi/nbdt
    interp_sp = Lagrange(11)
    interp_v = Lagrange(11)
    
    adv = Advection(
               (mesh_sp, mesh_v),
               [interp_sp, interp_v],
               dt,
               [([1, 2], 1, 1, true), ([2, 1], 1, 2, true)];
               tab_coef = magicsplit(dt),
           )
    
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

end

run_rotation()

```

## Vlasov/poisson

For Vlasov Poisson equation 

```@example landau

using SemiLagrangian
using Plots

function run_simulation(nbdt, sz, dt, interp, tab_coef)
    
    epsilon = 0.001

    xmin, xmax, nx = 0., 4π, sz[1]
    vmin, vmax, nv = -6., 6., sz[2]

    mesh_x = UniformMesh(xmin, xmax, nx)
    mesh_v = UniformMesh(vmin, vmax, nv)

    states = [([1, 2], 1, 1, true), ([2, 1], 1, 2, true)]

    adv = Advection((mesh_x, mesh_v), [interp, interp], dt, states; 
        tab_coef, timeopt = NoTimeOpt)
    
    kx = 0.5 
    fct_x(x) = epsilon * cos(kx * x) + 1
    fct_v(v) = exp(-v^2 / 2) / sqrt(2π)

    lgn_x = fct_x.(mesh_x.points)
    lgn_v = fct_v.(mesh_v.points)

    data = dotprod((lgn_x, lgn_v))

    pvar = getpoissonvar(adv)

    advd = AdvectionData(adv, data, pvar)

    time = Float64[]
    el = Float64[]
    for i = 1:nbdt
        while advection!(advd) end
        push!(time, advd.time_cur)
        push!(el, compute_ee(advd))
    end
    return time, el
end

nbdt = 1000
sz = (64, 64)
dt = 0.1
interp = Lagrange(9, Float64)
tab_coef = strangsplit(dt)
time, el = run_simulation( nbdt, sz, dt, interp, tab_coef)
plot(time, 0.5 .* log.(el.^2))

```

