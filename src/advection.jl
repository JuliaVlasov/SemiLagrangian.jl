"""
    advection!(data, mesh, velocity, dt, interpolant; dim=1)

Advect `data` in place along dimension `dim` for one time step `dt`.
`mesh` describes the advected dimension. `velocity` may be a scalar, an array
with one value for each line perpendicular to `dim`, or an array with the same
shape as `data` for a pointwise velocity field.

The characteristic shift is `-dt * velocity / step(mesh)`. Call this function
repeatedly to compose directional advections or construct a time-splitting
scheme. Pointwise velocity fields require an `AbstractInterpolation`; other
interpolants support scalar and per-line velocities.
"""
function advection!(
    data::AbstractArray,
    mesh::UniformMesh,
    velocity,
    dt::Real,
    interpolant;
    dim::Integer = 1,
)
    ndims(data) > 0 || throw(ArgumentError("data must have at least one dimension"))
    1 <= dim <= ndims(data) ||
        throw(ArgumentError("dim=$dim must be between 1 and $(ndims(data))"))
    size(data, dim) == length(mesh) ||
        throw(DimensionMismatch("size(data, $dim) must equal length(mesh)"))

    transverse_axes = Tuple(i for i in 1:ndims(data) if i != dim)
    transverse_shape = Tuple(size(data, i) for i in transverse_axes)
    transverse_indices = CartesianIndices(transverse_shape)
    npoints = size(data, dim)
    input = similar(data, eltype(data), (npoints,))
    output = similar(input)

    velocity_mode = if velocity isa Real
        :scalar
    elseif velocity isa AbstractArray{<:Real}
        if size(velocity) == size(data)
            :pointwise
        elseif size(velocity) == transverse_shape
            :per_line
        else
            throw(DimensionMismatch(
                "velocity must be a scalar, have transverse shape $transverse_shape, " *
                "or match data shape $(size(data)); got $(size(velocity))",
            ))
        end
    else
        throw(ArgumentError("velocity must be a real scalar or real array"))
    end

    per_line_velocity = velocity_mode == :per_line ? vec(velocity) : nothing
    scale = -dt / step(mesh)

    for (line, index) in enumerate(transverse_indices)
        transverse_index = Tuple(index)
        selectors = ntuple(
            i -> i == dim ? Colon() : transverse_index[i < dim ? i : i - 1],
            ndims(data),
        )
        data_line = view(data, selectors...)
        copyto!(input, data_line)
        if velocity_mode == :scalar
            _interpolate_advection_line!(output, input, scale * velocity, interpolant)
        elseif velocity_mode == :per_line
            _interpolate_advection_line!(
                output,
                input,
                scale * per_line_velocity[line],
                interpolant,
            )
        else
            shifts = scale .* view(velocity, selectors...)
            _interpolate_variable_advection_line!(output, input, shifts, interpolant)
        end
        copyto!(data_line, output)
    end

    return data
end

"""
    AdvectionOperator(mesh, velocity, interpolant; dim=1)

Describe one directional advection operator for use with [`split_advection!`](@ref).
"""
struct AdvectionOperator{M<:UniformMesh,V,I}
    mesh::M
    velocity::V
    interpolant::I
    dim::Int
end

function AdvectionOperator(mesh::UniformMesh, velocity, interpolant; dim::Integer = 1)
    return AdvectionOperator(mesh, velocity, interpolant, Int(dim))
end

"""
    split_advection!(data, first, second, coefficients)

Apply a two-operator splitting in place. Stages alternate between `first` and
`second`, beginning with `first`; `coefficients` are the time steps for the
individual stages, such as the result of `strangsplit(dt)`.
"""
function split_advection!(
    data::AbstractArray,
    first::AdvectionOperator,
    second::AdvectionOperator,
    coefficients::AbstractVector{<:Real},
)
    isempty(coefficients) && throw(ArgumentError("coefficients must not be empty"))

    for (stage, coefficient) in enumerate(coefficients)
        operator = isodd(stage) ? first : second
        advection!(
            data,
            operator.mesh,
            operator.velocity,
            coefficient,
            operator.interpolant;
            dim = operator.dim,
        )
    end

    return data
end

_interpolate_advection_line!(output, input, shift, interpolant::AbstractInterpolation) =
    interpolate!(output, input, shift, interpolant)

_interpolate_advection_line!(output, input, shift, interpolant) =
    interpolate!(output, interpolant, input, shift)

function _interpolate_variable_advection_line!(
    output,
    input,
    shifts,
    interpolant::AbstractInterpolation,
)
    return interpolate!(output, input, i -> (shifts[i],), [interpolant])
end

function _interpolate_variable_advection_line!(
    output,
    input,
    shifts,
    interpolant::AbstractShiftInterpolation,
)
    throw(ArgumentError(
        "pointwise velocity fields are not supported by this interpolant; " *
        "use a scalar or per-line velocity",
    ))
end

function _interpolate_variable_advection_line!(output, input, shifts, interpolant)
    throw(ArgumentError(
        "pointwise velocity fields require an AbstractInterpolation; " *
        "use a scalar or per-line velocity with this interpolant",
    ))
end
