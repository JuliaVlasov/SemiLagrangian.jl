using FFTW, LinearAlgebra

export Bspline

struct Bspline <: InterpolationType

    p::Int

    function Bspline(p::Int)

        @assert (p & 1 == 1)
        new(p)

    end
end


"""
    bspline(p, j, x)

Return the value at x in [0,1[ of the B-spline with
integer nodes of degree p with support starting at j.
Implemented recursively using the de Boor's recursion formula
using the [De Boor's Algorithm](https://en.wikipedia.org/wiki/De_Boor%27s_algorithm)

```math
B_{i,0}(x) := \\left\\{
\\begin{matrix}
1 & \\mathrm{if}  \\quad t_i ≤ x < t_{i+1} \\\\
0 & \\mathrm{otherwise}
\\end{matrix}
\\right.
```

```math
B_{i,p}(x) := \\frac{x - t_i}{t_{i+p} - t_i} B_{i,p-1}(x)
+ \\frac{t_{i+p+1} - x}{t_{i+p+1} - t_{i+1}} B_{i+1,p-1}(x).
```

"""
function bspline(p::Int, j::Int, x::Float64)

    if p == 0
        if j == 0
            return 1.0
        else
            return 0.0
        end
    else
        w::Float64 = (x - j) / p
        w1::Float64 = (x - j - 1) / p
    end
    (w * bspline(p - 1, j, x) + (1 - w1) * bspline(p - 1, j + 1, x))

end

export PeriodicAdvection

"""
    advection! = PeriodicAdvection( mesh, bspl )

Type to perform 1d advection on periodic domain. Bspline interpolation
is used combined with fft transform to solve the system. Bspline degree
must be odd.

```@example

p = 5
nx = 128

x1min, x1max = -10, 10

mesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)

advection! = PeriodicAdvection( mesh1, Bspline(p) )
```

advection! computes the interpolating spline of degree p of odd
degree of first dimension of array f on a periodic uniform mesh, at
all points x-alpha. f type is Array{Float64,2}.

"""
struct PeriodicAdvection <: AbstractAdvection

    p::Int
    mesh::UniformMesh
    modes::Vector{ComplexF64}
    eig_bspl::Vector{ComplexF64}
    eigalpha::Vector{ComplexF64}

    function PeriodicAdvection(mesh::UniformMesh, bspl::Bspline)

        n = mesh.length
        modes = zeros(Float64, n)
        modes .= 2π .* (0:n-1) ./ n
        eig_bspl = zeros(Complex{Float64}, n)
        eigalpha = zeros(Complex{Float64}, n)
        eig_bspl .= bspline(bspl.p, -div(bspl.p + 1, 2), 0.0)
        for j = 1:div(bspl.p + 1, 2)-1
            eig_bspl .+= (bspline(bspl.p, j - div(bspl.p + 1, 2), 0.0) .* 2 .*
                          cos.(j .* modes))
        end

        new(bspl.p, mesh, modes, eig_bspl, eigalpha)

    end

end

function (adv::PeriodicAdvection)(
    f::Array{Float64,2},
    v::Vector{Float64},
    dt::Float64,
)
    nv = length(v)
    p = adv.p
    delta = adv.mesh.step

    @assert size(f)[2] == nv

    f̂ = fft(f, 1)

    for j = 1:nv
        alpha = v[j] * dt / delta
        ishift = floor(-alpha)
        beta = -ishift - alpha
        fill!(adv.eigalpha, 0.0im)
        for i = -div(p - 1, 2):div(p + 1, 2)
            adv.eigalpha .+= (bspline(p, i - div(p + 1, 2), beta) .*
                              exp.((ishift + i) * 1im .* adv.modes))
        end

        f̂[:, j] .= f̂[:, j] .* adv.eigalpha ./ adv.eig_bspl

    end

    f .= real(ifft(f̂, 1))

end


function (adv::PeriodicAdvection)(
    f::Array{ComplexF64,2},
    v::Vector{Float64},
    dt::Float64,
)

    nv = length(v)
    dx = adv.mesh.step

    @assert size(f)[2] == nv

    fft!(f, 1)

    @inbounds for j = 1:nv
        alpha = dt * v[j] / dx
      # compute eigenvalues of cubic splines evaluated at displaced points
        ishift = floor(-alpha)
        beta = -ishift - alpha
        fill!(adv.eigalpha, 0.0im)
        for i = -div(adv.p - 1, 2):div(adv.p + 1, 2)
            adv.eigalpha .+= (bspline(adv.p, i - div(adv.p + 1, 2), beta) .*
                              exp.((ishift + i) * 1im .* adv.modes))
        end

      # compute interpolating spline using fft and
      # properties of circulant matrices

        f[:, j] .*= adv.eigalpha ./ adv.eig_bspl

    end

    ifft!(f, 1)

end
