# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Julia 1.8.5
#     language: julia
#     name: julia-1.8
# ---

# # Semi-Lagrangian inverse problem

using PyPlot
using FFTW
using LinearAlgebra

# +
struct Mesh
    
    nx :: Int
    nv :: Int
    lx :: Float64
    lv :: Float64
    dx :: Float64
    dv :: Float64
    x :: Vector{Float64}
    v :: Vector{Float64}
    
    function Mesh(lx, nx, lv, nv)
        x = LinRange(0, lx, nx+1)[1:end-1] 
        dx = lx / nx
        v = LinRange(-lv, lv, nv+1)[1:end-1] 
        dv = lv / nv
        new(nx, nv, lx, lv, dx, dv, x, v)
    end

end
# -

function bspline(p::Int, j::Int, x::Float64)
   if p == 0
       if j == 0
           return 1.0
       else
           return 0.0
       end
   else
       w = (x - j) / p
       w1 = (x - j - 1) / p
   end
   return (w * bspline(p - 1, j, x) + (1 - w1) * bspline(p - 1, j + 1, x))
end

# +
function advection!(f::Matrix{Float64}, v::Vector{Float64}, delta, dt::Float64)

   p = 3
   n = size(f, 1)
   modes = [2π * i / n for i in 0:n-1]
   eig_bspl = zeros(Float64, n)
   eig_bspl .= bspline(p, -div(p+1,2), 0.0)
   for i in 1:div(p+1,2)-1
      eig_bspl .+= bspline(p, i - (p+1)÷2, 0.0) * 2 .* cos.(i * modes)
   end
   eigalpha = zeros(ComplexF64, n)
   f̂ = zeros(ComplexF64, size(f))
   f̂ .= f
   fft!(f̂,1)
    
   for j in eachindex(v)
      @inbounds alpha = dt * v[j] / delta
      ishift = floor(-alpha)
      beta   = -ishift - alpha
      fill!(eigalpha,0.0im)
      for i in -div(p-1,2):div(p+1,2)
         eigalpha .+= (bspline(p, i-div(p+1,2), beta) 
                        .* exp.((ishift+i) * 1im .* modes))
      end      
      @inbounds f̂[:,j] .*= eigalpha ./ eig_bspl
        
   end

   ifft!(f̂,1)
   f .= real(f̂)
    
end            
# -

compute_rho(mesh::Mesh, f) = vec(1.0 .- mesh.dv * sum(f, dims=2))


function compute_e_from_rho(mesh::Mesh, ρ)
    nx = mesh.nx
    ρ̂ = fft(ρ)
    ê = 1im .* zero(ρ̂)
    ik = 1im .* 2π ./ mesh.lx .* collect(fftfreq(nx, nx))
    ê[2:end] .= - ρ̂[2:end] ./ ik[2:end]
    ê[1] = 0.0
    return real(ifft(ê))
end


# +
compute_e(mesh::Mesh, f) = compute_e_from_rho(mesh, compute_rho(mesh, f))

electric_energy(mesh::Mesh, e) = 0.5 * sum( e.^2 ) * mesh.dx
# -


nx, nv = 2*163, 543
lx, lv = 20π/3, 12.0
mesh = Mesh(lx, nx, lv, nv)
num_steps = 300
deltat  = 0.1
fv = (0.9 .* exp.(-mesh.v.^2)./sqrt(π)  .+ 0.1 .* exp.(-2 .* (mesh.v .- 4.5).^2) ./ sqrt(π/2))
f_eq = zeros(nx, nv)
f_eq = fv' .* ones(nx)
imshow(f_eq)

h = zero(mesh.x)
f0 = zeros(nx, nv)
f0 .= (1.0 .+ 0.03 .* cos.(0.3 .* mesh.x)) .* f_eq 
surf(f0)

# Compute $ρ(x)$

ρ = compute_rho(mesh, f0)

plot(mesh.x, ρ);

e = compute_e_from_rho(mesh, ρ)

plot(mesh.x, e)

# ## Forward problem

electric_energy(mesh, e)

# +
using ProgressMeter

function run_forward(mesh, f_eq, h)
    f = zeros(mesh.nx, mesh.nv)
    f .= (1.0 .+ 0.03 .* cos.(0.3 .* mesh.x)) .* f_eq
    t = 0.0
    es = Float64[]
    ts = Float64[]
    fstar = Matrix{Float64}[]
    estar = Vector{Float64}[]
    ft = zeros(mesh.nv, mesh.nx)
    e = zeros(mesh.nx)
    e_total = zeros(mesh.nx)
    
    @showprogress 1 for i in 1:num_steps
        e .= compute_e(mesh, f)
        ee = electric_energy(mesh, e)
        push!(es, ee)
        push!(ts, t)
        advection!(f, mesh.v, mesh.dx, 0.5deltat)
        push!(fstar, copy(f))
        e_total .= compute_e(mesh, f) .+ h   
        push!(estar, copy(e_total))
        transpose!(ft, f)
        advection!(ft, e_total, mesh.dv, deltat)
        transpose!(f, ft)
        advection!(f, mesh.v, mesh.dx, 0.5deltat)
        t += deltat
    end
    return ts, es, fstar, estar
end
# -


@time ts, es, fstar, estar = run_forward(mesh, f_eq, h);

plot(ts, log.(es))

# ## Backward problem

function shift_by_n(mesh::Mesh, dt, g, eh)
    out = zero(g)
    for i in 1:mesh.nx
        n = floor(Int, -dt*eh[i]/mesh.dv)
        out[i, :] .= circshift(g[i, :], -n)
    end
    return out
end

# +
function run_backward(mesh::Mesh, f, fstar, estar)
    ft = zeros(mesh.nv, mesh.nx)
    @showprogress 1 for i in reverse(1:num_steps)
        advection!(f, mesh.v, mesh.dx, -0.5deltat)
        transpose!(ft, f)
        advection!(ft, estar[i], mesh.dv, -deltat)
        transpose!(f, ft)
        advection!(f, mesh.v, mesh.dx, -0.5deltat)
    end

    return f
        
end


# +
gt = 2 .* (fstar[end] - f_eq) / deltat

@time gs = run_backward(mesh, gt, fstar, estar)
# -

energy(mesh, f, f_eq) = sum((f .- f_eq).^2) * mesh.dx * mesh.dv

ee = [electric_energy(mesh, e) for e in estar]
hist_j = [energy(mesh, f, f_eq) for f in fstar]
plot_all(mesh, gs, ee, h, hist_j)

# +
# gradient
f_ss_s = [shift_by_n(0 * deltat, f, eh) for (f, eh) in zip(f_stars,e_stars)]
f_ss_d = [(circshift(f, (0,1)) - circshift(f, (0,-1)))/ 2 / hv for f in f_ss_s]
gradient = sum([sum(f, dims=2) * hv for (f, g) in zip(f_ss_s, g_starstars)], dims=1) * deltat^2

plot(xs, gradient, label="gradient")
legend()
# -

plot_all(fs[end], ee, h, hist_j)

plot_all(g_starstars[end], ee, h, hist_j)



# +
function plot_all(mesh::Mesh, f, ee, h, hist_j)
    
    imshow(f, extent=[mesh.x[begin], mesh.x[end], mesh.v[begin], mesh.v[end]])
    colorbar()
    xlabel("x")
    ylabel("v")

    subplot(2,2,2)
    plot(mesh.x, compute_rho(mesh, f), label="rho")
    plot(mesh.x, compute_e(mesh, f), label="e")
    plot(mesh.x, h, label="H")
    xlabel("x")
    ylabel("E")
    legend()

    subplot(2,2,3)
    semilogy(1:length(ee) .* deltat, ee)
    xlabel("t")
    ylabel("electric energy")
    
    subplot(2,2,4)
    semilogy(1:length(hist_J), hist_J)
    xlabel("iterations")
    ylabel("J")
    
end

j(mesh, f) = sum((f[end] .- f_eq).^2) * mesh.dx * mesh.dv

function gradient_fd(h, d)
    eps = 1e-10
    f0, _, _, _ = run_forward(h)
    f1, _, _, _ = run_forward(h + eps*d)
    return (j(f1) - j(f0))/eps
end
