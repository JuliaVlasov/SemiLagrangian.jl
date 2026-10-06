struct Polynomial{T}
    coefficients::Vector{T}
    function Polynomial{T}(coefficients::AbstractVector) where {T}
        isempty(coefficients) && return new{T}(T[zero(T)])
        values = T.(coefficients)
        while length(values) > 1 && iszero(last(values))
            pop!(values)
        end
        return new{T}(values)
    end
end

Polynomial(coefficients::AbstractVector{T}) where {T} = Polynomial{T}(coefficients)
Base.eltype(::Type{Polynomial{T}}) where {T} = T
Base.convert(::Type{Polynomial{T}}, polynomial::Polynomial) where {T} =
    Polynomial{T}(polynomial.coefficients)
Base.zero(::Type{Polynomial{T}}) where {T} = Polynomial{T}(T[zero(T)])
Base.one(::Type{Polynomial{T}}) where {T} = Polynomial{T}(T[one(T)])
Base.zero(polynomial::Polynomial) = zero(typeof(polynomial))
Base.one(polynomial::Polynomial) = one(typeof(polynomial))
Base.iszero(polynomial::Polynomial) = all(iszero, polynomial.coefficients)
Base.copy(polynomial::Polynomial) = Polynomial(copy(polynomial.coefficients))
Base.length(polynomial::Polynomial) = length(polynomial.coefficients)
Base.:(==)(a::Polynomial, b::Polynomial) = a.coefficients == b.coefficients
degree(polynomial::Polynomial) = iszero(polynomial) ? -1 : length(polynomial) - 1

function (polynomial::Polynomial{T})(x::Number) where {T}
    return evalpoly(x, polynomial.coefficients)
end

function (polynomial::Polynomial)(x::Polynomial)
    result = zero(x)
    for coefficient in Iterators.reverse(polynomial.coefficients)
        result = result * x + coefficient
    end
    return result
end

function Base.:+(a::Polynomial, b::Polynomial)
    T = promote_type(eltype(a.coefficients), eltype(b.coefficients))
    result = zeros(T, max(length(a), length(b)))
    @inbounds for i in eachindex(a.coefficients)
        result[i] += a.coefficients[i]
    end
    @inbounds for i in eachindex(b.coefficients)
        result[i] += b.coefficients[i]
    end
    return Polynomial(result)
end

Base.:-(a::Polynomial) = Polynomial(-a.coefficients)
Base.:-(a::Polynomial, b::Polynomial) = a + (-b)

function Base.:+(polynomial::Polynomial, scalar::Number)
    return polynomial + Polynomial([scalar])
end
Base.:+(scalar::Number, polynomial::Polynomial) = polynomial + scalar
Base.:-(polynomial::Polynomial, scalar::Number) = polynomial + (-scalar)
Base.:-(scalar::Number, polynomial::Polynomial) = (-polynomial) + scalar

function Base.:*(a::Polynomial, b::Polynomial)
    T = promote_type(eltype(a.coefficients), eltype(b.coefficients))
    result = zeros(T, length(a) + length(b) - 1)
    @inbounds for i in eachindex(a.coefficients), j in eachindex(b.coefficients)
        result[i + j - 1] += a.coefficients[i] * b.coefficients[j]
    end
    return Polynomial(result)
end

Base.:*(polynomial::Polynomial, scalar::Number) =
    Polynomial(polynomial.coefficients .* scalar)
Base.:*(scalar::Number, polynomial::Polynomial) = polynomial * scalar
Base.:/(polynomial::Polynomial, scalar::Number) =
    Polynomial(polynomial.coefficients ./ scalar)
function Base.:^(polynomial::Polynomial, n::Integer)
    n >= 0 || throw(DomainError(n, "polynomial exponent must be nonnegative"))
    result = one(polynomial)
    factor = polynomial
    exponent = n
    while exponent > 0
        isodd(exponent) && (result *= factor)
        exponent >>= 1
        exponent > 0 && (factor *= factor)
    end
    return result
end

function derivative(polynomial::Polynomial)
    length(polynomial) <= 1 && return zero(typeof(polynomial))
    return Polynomial([
        i * polynomial.coefficients[i + 1] for i in 1:(length(polynomial) - 1)
    ])
end
