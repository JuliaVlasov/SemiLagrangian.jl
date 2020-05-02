"""
    advection! = PeriodicAdvection( mesh, lagrange )

Type to perform 1d advection on periodic domain. Bspline or Lagrange 
interpolation is used.

```@example

p = 5
n1, n2 = 128, 128

x1min, x1max = -10, 10
x2min, x2max = -10, 10

mesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)
mesh2 = UniformMesh(x2min, x2max, n2; endpoint=false)

advection1! = PeriodicAdvection( mesh1, BSpline(p) )
advection2! = PeriodicAdvection( mesh2, Lagrange(p) )
```

advection1! computes the interpolating spline of degree p of odd
degree of first dimension of array f on a periodic uniform mesh, at
all points x-alpha. f type is Array{ComplexF64,2}.

advection2! computes the Lagrange interpolation with stencil p
on first dimension of array f on a periodic uniform mesh, at
all points x-alpha. f type is Array{Float64,2}.

"""
struct LagrangePeriodicAdvection <: AbstractAdvection

    mesh::UniformMesh
    lag::Lagrange
    fp::Vector{Float64}

    function LagrangePeriodicAdvection(mesh::UniformMesh, lag::Lagrange)

        n = mesh.length

        new(mesh, lag, zeros(n))

    end

end


function (adv::LagrangePeriodicAdvection)(
    f::Array{Float64,2},
    v::Vector{Float64},
    dt::Float64,
)
    nv = length(v)
    stencil = adv.lag.stencil
    delta = adv.mesh.step

    @assert size(f)[2] == nv

    for j in eachindex(v)
        alpha = - v[j] * dt / delta

        fi = view(f,:,j)
        lagrange_interpolation_1d_disp_fixed_periodic(fi, adv.fp, alpha, stencil)
        f[:,j] .= adv.fp

    end

end


