# Rotation

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

```@autodocs
Modules = [SemiLagrangian]
Order   = [:type, :function]
Pages   = ["rotation.jl"]
```
