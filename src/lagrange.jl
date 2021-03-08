

"""
    _getpolylagrange(k::Int64, order::Int64, origin::Int64)

Function that return the k-th Lagrange Polynomial of a certain order. Coefficients are rational then the return is exact. The polynomial is equal to :
``\\prod_{i=0,\\ i \\neq k}^{order} \\frac{x - i - origin}{k - i}``

# Arguments
- `k::Int64` : number of the Polynomial, `k` must be between `0` and `order` (`0<= k <= order`).
- `order::Int64` : order of the polynomial.
- `origin::Int64` : origin of the first indice.

# Returns
- `Polynomial{Rational{BigInt}}` : the k-th Lagrange polynomial of order `order`

# Throws
- `DommaineError` : when `0 <= k <= order` is `false` or when N ∉ {BInt64, BigInt}
"""
function _getpolylagrange(k::Int, order::Int, origin::Int) where {N<:Integer}
    0 <= k <= order || throw(DomainError("the constant 0 <= k <= order is false"))
    # the computed is made with big rational
    result = Polynomials.Polynomial([big(1//1)])
    for l=0:order
        if l != k
            result *= Polynomials.Polynomial([-(l+origin),1//1] .// (k-l) )
        end
    end
    return result
end

"""
    Lagrange{T, edge, order, N} <: AbstractInterpolation{T, edge, order}
    Lagrange(order, T::DataType=Float64; edge::EdgeType=CircEdge)

Type containing Lagrange Polynomials coefficients for Lagrange interpolation

# Type parameters
- `T` : the type of data that is interpolate
- `edge::EdgeType` : type of edge traitment
- `order::Int`: order of lagrange interpolation

# Implementation :
- `tabfct::Vector{Polynomial{T}}` : vector of all lagrange polynomial, per example the k-th Lagrange polynomial for the designed order is tabfct[k+1]

# Arguments : 
- `order::Int` : the order of interpolation
- `[T::DataType=Float64]` : The type values to interpolate 

# Keywords arguments :
- `edge::EdgeType=CircEdge` : type of edge traitment

"""
struct Lagrange{T, edge, order} <: AbstractInterpolation{T, edge, order}
    tabfct::Vector{Polynomial{T}}
    function Lagrange(order::Int, T::DataType=Float64; edge::EdgeType=CircEdge) 
        origin = -div(order,2)
        tabfct_rat = collect([_getpolylagrange( i, order, origin) for i=0:order])

        new{T, edge, order}(convert.(Polynomial{T}, tabfct_rat)) 
    end
end


struct Lagrange2d{T, edge, order} <: AbstractInterpolation2d{T, edge, order}
    l1d::Lagrange{T, edge, order}
    function Lagrange2d(order, T::DataType=Float64; edge::EdgeType=CircEdge)
        new{T, edge, order}(Lagrange(order, T, edge=edge))
    end
end

gettabfct(interp::Lagrange2d)=interp.l1d.tabfct
