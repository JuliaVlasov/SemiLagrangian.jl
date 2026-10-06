using FastInterpolations
using FFTW
using LinearAlgebra
using Plots
using DispersionRelations
import Statistics: mean
using .Threads

struct UniformMesh
    xmin::Float64
    xmax::Float64
    nx::Int
    dx::Float64
    x::Vector{Float64}
    function UniformMesh(xmin, xmax, nx)
        dx = (xmax - xmin) / nx
        x = LinRange(xmin, xmax, nx + 1)
        return new(xmin, xmax, nx, dx, x)
    end
end

function compute_rho!(
        rho::Vector{Float64},
        meshv::UniformMesh,
        f::Array{Float64, 2}
    )

    dv = meshv.dx
    rho .= dv .* vec(sum(f, dims = 1))
    rho .-= mean(rho)
    rho[nx+1] = rho[1]
end

struct PoissonSolver

    modes :: Vector{Float64}
    rhok :: Vector{ComplexF64}

    function PoissonSolver(meshx)

        nx = meshx.nx
        modes = 2π / (meshx.xmax - meshx.xmin) * collect(fftfreq(nx, nx))
        modes[1] = 1.0
        rhok = zeros(ComplexF64, nx)
        new( modes, rhok )

    end

end

function compute_e!(e::Vector{Float64}, solver::PoissonSolver, rho::Vector{Float64})
    solver.rhok .= -1im .* fft(rho[1:end-1]) ./ solver.modes 
    e[1:nx] .= real(ifft(solver.rhok))
    e[nx+1] = e[1]
end

function advection_x!(
        f::Matrix{Float64},
        mesh::UniformMesh, v::Vector{Float64},
        dt::Float64
    )

    nx, dx = mesh.nx, mesh.dx
    xi = 1.0:nx+1
    cache = CubicSplineCache(1.0:nx+1; bc = PeriodicBC())

    @threads for j in eachindex(v)
        alpha = - dt * v[j] / dx
        fi = view(f, :, j)
        xp = xi .+ alpha
        fi[nx+1] = fi[1]
        f[:, j] .= cubic_interp(cache, fi, xp) 
    end

end

function advection_v!(
        f::Matrix{Float64},
        mesh::UniformMesh, v::Vector{Float64},
        dt::Float64
    )

    nx, dx = mesh.nx, mesh.dx
    xi = 1.0:nx+1
    cache = CubicSplineCache(1.0:nx+1; bc=ZeroSlopeBC())  

    @threads for j in eachindex(v)
        alpha = - dt * v[j] / dx
        fi = view(f, :, j)
        xp = xi .+ alpha
        xp = max.(xp, 1)
        xp = min.(xp, nx+1)
        fi[nx+1] = fi[1]
        f[:, j] .= cubic_interp(cache, fi, xp) 
    end
end

function landau_fast(nx, nv, dt, nt::Int64)

    # Set grid
    eps, kx = 0.001, 0.4
    xmin, xmax = 0.0, 2π / kx
    vmin, vmax = -6.0, 6.0
    meshx = UniformMesh(xmin, xmax, nx)
    meshv = UniformMesh(vmin, vmax, nv)

    x = meshx.x
    v = meshv.x
    dx = meshx.dx

    f = zeros(Float64, (nx+1, nv+1))
    f .= (1.0 .+ eps * cos.(kx * x)) / sqrt(2π) .* transpose(exp.(-0.5 * v .^ 2))
    fᵗ = zeros(Float64, (nv+1, nx+1))
    permutedims!(fᵗ, f, [2,1])

    poisson = PoissonSolver(meshx)

    rho = zeros(nx+1)
    e = zeros(nx+1)

    compute_rho!(rho, meshv, fᵗ)
    compute_e!(e, poisson, rho)
    ℰ = Float64[sum(e .* e) * dx]
    t = Float64[0.0]

    for it in 1:nt
        advection_x!(f, meshx, v, 0.5dt)
        permutedims!(fᵗ, f, [2,1])
        compute_rho!(rho, meshv, fᵗ)
        compute_e!(e, poisson, rho)
        push!(ℰ, sum(e .* e) * dx )
        push!(t, (it - 0.5) * dt)
        advection_v!(fᵗ, meshv, e, dt)
        permutedims!(f, fᵗ, [2,1])
        advection_x!(f, meshx, v, 0.5dt)
    end

    return t, ℰ

end

nx, nv = 256, 512
dt, nt = 0.1, 1000
t, nrj = landau_fast(nx, nv, dt, 1)
landau_fast(nx, nv, dt, 1) # warmup
@time t, nrj = landau_fast(nx, nv, dt, nt)
plot(t, nrj; label = "|E|²", yaxis = :log)
line, ω, = fit_complex_frequency(t, nrj)
plot!(t, line, yaxis = :log, label = "$(imag(ω/2))")
title!("α = 0.001, k = 0.4")

# about fit_complex_frequency
# if E = |exp.(-im * (a + im * b) * t)|^2
# then ω = 2 * (a + im * b)
