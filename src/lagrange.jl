

"""
    _getpolylagrange(k::Int64, order::Int64, origin::Int64, N::DataType)

Function that return the k-th Lagrange Polynomial of a certain order. Coefficients are rational then the return is exact. The polynomial is equal to :
``\\prod_{i=0,\\ i \\neq k}^{order} \\frac{x - i - origin}{k - i}``

# Arguments
- `k::Int64` : number of the Polynomial, `k` must be between `0` and `order` (`0<= k <= order`).
- `order::Int64` : order of the polynomial.
- `origin::Int64` : origin of the first indice.
- `N::DataType` : type of the Integer that is the base of the rational type, in fact Int64 or BigInt. BigInt is needed for order greater or equal to 21.

# Returns
- `Polynomial{Rational{N}}` : the k-th Lagrange polynomial of order `order`

# Throws
- `DommaineError` : when `0 <= k <= order` is `false` or when N ∉ {BInt64, BigInt}
"""
function _getpolylagrange(k::Int64, order::Int64, origin::Int64, fact::N) where {N}
    0 <= k <= order || throw(DomainError("the constant 0 <= k <= order is false"))
    N <: Union{BigInt,Int64} || throw(DomainError(N, "N must be Int64 or BigInt"))
    # the computed is made with big integer
    result = Polynomials.Polynomial([big(fact//1)])
    for l=0:order
        if l != k
            result *= Polynomials.Polynomial([-(l+origin),1//1] .// (k-l) )
         end
    end
    return Polynomials.Polynomial(N.(result.coeffs))
end

"""
    Lagrange{T, iscirc, order, N} <: InterpolationType{T, iscirc, order}
    Lagrange(order, T::DataType=Float64; iscirc::Bool=true)

Get Lagrange Polynomials coefficients for Lagrange interpolation

# Type parameters
- `T` : the type of data that is interpolate
- `iscirc::Bool` : true if function is circular
- `order::Int`: order of lagrange interpolation
- `N` : type of integer, in fact Int64 or BigInt that is used to store lagrange polynomial

# Implementation :
- `fact_order::N` : factorial of the order
- `lagpol:Vector{Polynomial{N}}` : vector of all lagrange polynomial, per example the k-th Lagrange polynomial for the designed order is lagpol[k+1]/fact_order
"""
struct Lagrange{T, iscirc, order, N} <: InterpolationType{T, iscirc, order}
    fact_order::N
    lagpol::Vector{Polynomial{N}}
    function Lagrange(order, T::DataType=Float64; iscirc::Bool=true) 
        type = order <= 20 ? Int64 : BigInt
        fact_order = factorial(type(order))
        origin = -div(order,2)
        lagpol = collect([_getpolylagrange( i, order, origin, fact_order) for i=0:order])
        new{T, iscirc, order, type}(fact_order, lagpol) 
    end
end
@inline get_order(lag::Lagrange{T,iscirc, order}) where{T, iscirc, order}= order
@inline get_type(lag::Lagrange{T, isc, order, N}) where{T, isc, order, N}="Lagrange{$T, $isc, $order, $N}"
@inline get_precal(lag::Lagrange{T}, decf) where{T}=@inbounds [T(fct(decf))/lag.fact_order for fct in lag.lagpol]
@inline get_precal!(v::Vector{T}, lag::Lagrange{T},decf) where{T}=@inbounds v .= get_precal(lag, decf)
@inline sol(lag::Lagrange, b)=b
@inline isbspline(_::Lagrange)=false
Base.show(io::IO, lag::Lagrange)=print(io, get_type(lag))
