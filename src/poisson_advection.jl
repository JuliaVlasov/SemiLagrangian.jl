"""
$(TYPEDEF)

    PoissonAdvectionOperator(meshx, meshv, interpolant = Spectral(length(meshv)))

Self-consistent velocity advection for periodic 1D-1V Vlasov-Poisson simulations.
Data is a `Float64` matrix with position along dimension 1 and velocity along
dimension 2. Both meshes must have positive spacing.

Each stage integrates the distribution over velocity, subtracts the mean
density, solves `dE/dx = rho` with zero mean electric field, and advects in
velocity using `E`. Density, field, FFT plans, and line buffers are reused.
The default spectral interpolation also uses cached FFT plans; other
interpolants may allocate internally. This mutable workspace is not thread-safe.

The latest density and field are available as `operator.rho` and `operator.field`.
Call [`update_field!`](@ref) to refresh them without advancing the distribution.
"""
struct PoissonAdvectionOperator{I, P, Q, L, R}
    meshx::UniformMesh{Float64}
    meshv::UniformMesh{Float64}
    interpolant::I
    rho::Vector{Float64}
    field::Vector{Float64}
    fourier::Vector{ComplexF64}
    inverse_modes::Vector{ComplexF64}
    forward::P
    backward::Q
    input::Vector{Float64}
    output::Vector{Float64}
    line_forward::L
    line_backward::R
end

function PoissonAdvectionOperator(
        meshx::UniformMesh{Float64},
        meshv::UniformMesh{Float64},
        interpolant = Spectral(length(meshv)),
    )
    nx, nv = length(meshx), length(meshv)
    nx > 1 && nv > 1 || throw(ArgumentError("both meshes need at least two points"))
    step(meshx) > 0 && step(meshv) > 0 ||
        throw(ArgumentError("mesh spacing must be positive"))
    if interpolant isa Union{Spectral, LagrangeFFT, CubicSpline}
        interpolant.nx == nv || throw(DimensionMismatch("interpolant length must equal velocity grid length"))
    end
    fourier = zeros(ComplexF64, nx)
    forward = plan_fft!(fourier)
    backward = plan_ifft!(fourier)
    modes = vec_k_fft(meshx)
    inverse_modes = ComplexF64[iszero(k) ? 0 : -im / k for k in modes]
    # A real field has no resolvable derivative at the even-grid Nyquist mode.
    iseven(nx) && (inverse_modes[nx ÷ 2 + 1] = 0)
    line_forward = interpolant isa Spectral ? plan_fft!(interpolant.ufft) : nothing
    line_backward = interpolant isa Spectral ? plan_ifft!(interpolant.ufft) : nothing
    return PoissonAdvectionOperator(
        meshx, meshv, interpolant, zeros(nx), zeros(nx), fourier,
        inverse_modes, forward, backward, zeros(nv), zeros(nv),
        line_forward, line_backward,
    )
end

"""
$(SIGNATURES)

Update the neutralized density and zero-mean electric field from `data`, without
advecting it. Reuses the operator's buffers and FFT plans; returns its field.
"""
function update_field!(operator::PoissonAdvectionOperator, data::AbstractMatrix{Float64})
    size(data) == (length(operator.meshx), length(operator.meshv)) ||
        throw(DimensionMismatch("data shape must match the position and velocity meshes"))
    fill!(operator.rho, 0)
    for j in axes(data, 2), i in axes(data, 1)
        operator.rho[i] += data[i, j]
    end
    operator.rho .*= step(operator.meshv)
    average = sum(operator.rho) / length(operator.rho)
    operator.rho .-= average
    operator.fourier .= operator.rho
    operator.forward * operator.fourier
    operator.fourier .*= operator.inverse_modes
    operator.backward * operator.fourier
    operator.field .= real.(operator.fourier)
    return operator.field
end

function _poisson_interpolate!(operator, interpolant, shift)
    return interpolate!(operator.output, operator.input, shift, interpolant)
end

function _poisson_interpolate!(operator, interpolant::Spectral, shift)
    n = length(operator.input)
    interpolant.ufft .= operator.input
    operator.line_forward * interpolant.ufft
    for j in 0:(n - 1)
        k = j <= n ÷ 2 ? j : j - n
        interpolant.ufft[j + 1] *= cis(2pi * k * shift / n)
    end
    operator.line_backward * interpolant.ufft
    operator.output .= real.(interpolant.ufft)
    return operator.output
end

"""
$(SIGNATURES)

Refresh the Poisson field from `data`, then advance the velocity direction by
`dt`. Returns `data`. The field is held fixed during this stage.
"""
function advection!(data::AbstractMatrix{Float64}, operator::PoissonAdvectionOperator, dt::Real)
    update_field!(operator, data)
    scale = -dt / step(operator.meshv)
    for i in axes(data, 1)
        for j in axes(data, 2)
            operator.input[j] = data[i, j]
        end
        _poisson_interpolate!(operator, operator.interpolant, scale * operator.field[i])
        for j in axes(data, 2)
            data[i, j] = operator.output[j]
        end
    end
    return data
end

function split_advection!(
        data::AbstractArray,
        first::Union{AdvectionOperator, PoissonAdvectionOperator},
        second::Union{AdvectionOperator, PoissonAdvectionOperator},
        coefficients::AbstractVector{<:Real},
    )
    isempty(coefficients) && throw(ArgumentError("coefficients must not be empty"))
    for (stage, coefficient) in enumerate(coefficients)
        advection!(data, isodd(stage) ? first : second, coefficient)
    end
    return data
end
