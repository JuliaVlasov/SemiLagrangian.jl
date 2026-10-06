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
        x = LinRange(xmin, xmax, nx + 1)[1:end-1]
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
    solver.rhok .= -1im .* fft(rho) ./ solver.modes 
    e .= real(ifft(solver.rhok))
end

function advection!(
        f::Matrix{Float64},
        mesh::UniformMesh, v::Vector{Float64},
        dt::Float64
    )

    nx, dx = mesh.nx, mesh.dx
    xi = 1.0:nx

    @threads for j in eachindex(v)
        alpha = - dt * v[j] / dx
        fi = view(f, :, j)
        xp = xi .+ alpha
        f[:, j] .= cubic_interp(xi, fi, xp; bc=PeriodicBC(endpoint=:exclusive))
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

    f = zeros(Float64, (nx, nv))
    f .= (1.0 .+ eps * cos.(kx * x)) / sqrt(2π) .* transpose(exp.(-0.5 * v .^ 2))
    fᵗ = zeros(Float64, (nv, nx))
    permutedims!(fᵗ, f, [2,1])

    poisson = PoissonSolver(meshx)

    rho = zeros(nx)
    e = zeros(nx)

    compute_rho!(rho, meshv, fᵗ)
    compute_e!(e, poisson, rho)
    ℰ = Float64[sum(e .* e) * dx]
    t = Float64[0.0]

    for it in 1:nt
        advection!(f, meshx, v, 0.5dt)
        permutedims!(fᵗ, f, [2,1])
        compute_rho!(rho, meshv, fᵗ)
        compute_e!(e, poisson, rho)
        push!(ℰ, sum(e .* e) * dx )
        push!(t, (it - 0.5) * dt)
        advection!(fᵗ, meshv, e, dt)
        permutedims!(f, fᵗ, [2,1])
        advection!(f, meshx, v, 0.5dt)
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
