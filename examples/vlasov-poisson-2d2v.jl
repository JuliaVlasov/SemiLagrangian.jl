using SemiLagrangian

mesh_sp = UniformMesh( 0.0, 4pi, 32)
mesh_v = UniformMesh( -6.0, 6.0, 32)

nbdt=50
dt = 0.1
epsilon = 0.5

interp = [ntuple(x -> Lagrange(7), 4)...]

tabst = [([3,4,1,2], 1,1, true), ([4,3,1,2], 1, 1, true),  ([1,2,4,3], 1, 2, true),([2,1,3,4], 1, 2, true)]

adv = Advection( (mesh_sp, mesh_sp, mesh_v, mesh_v), interp, dt, tabst)

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
println("$t\t$ee\t$ke\t$(ee+ke)")

for ind=1:nbdt
    while advection!(advdata) end
    t = Float32(dt*ind)
    compute_charge!(advdata)
    compute_elfield!(advdata)
    ee = compute_ee(advdata)
    ke = compute_ke(advdata)
    println("$t\t$ee\t$ke\t$(ee+ke)")
end
