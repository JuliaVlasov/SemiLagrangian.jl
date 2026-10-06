
function _getlagrangecoefficients(k::Int, order::Int, origin::Int)
    0 <= k <= order || throw(DomainError("the constant 0 <= k <= order is false"))
    coefficients = [big(1 // 1)]
    for l = 0:order
        if l != k
            denominator = big(k - l)
            constant = big(-(l + origin)) // denominator
            linear = big(1) // denominator
            product = zeros(Rational{BigInt}, length(coefficients) + 1)
            for i in eachindex(coefficients)
                product[i] += coefficients[i] * constant
                product[i + 1] += coefficients[i] * linear
            end
            coefficients = product
        end
    end
    return coefficients
end

struct LagrangePolynomial{T}
    coefficients::Vector{T}
end

(polynomial::LagrangePolynomial)(x) = evalpoly(x, polynomial.coefficients)

"""
$(TYPEDEF)

Type containing the Lagrange basis polynomials for interpolation.

# Type parameters
- `T` : the type of data that is interpolate
- `edge::EdgeType` : type of edge traitment
- `order::Int`: order of lagrange interpolation

# Implementation :
- `tabfct::Vector{LagrangePolynomial{T}}` : callable basis polynomials, with the k-th basis polynomial at `tabfct[k+1]`

# Arguments : 
- `order::Int` : the order of interpolation
- `[T::DataType=Float64]` : The type values to interpolate 

# Keywords arguments :
- `edge::EdgeType=CircEdge` : type of edge traitment

"""
struct Lagrange{T,edge,order} <: AbstractInterpolation{T,edge,order}

    tabfct::Vector{LagrangePolynomial{T}}

    function Lagrange(order::Int, T::DataType = Float64; edge::EdgeType = CircEdge)

        origin = -div(order, 2)

        tabfct_rat = [_getlagrangecoefficients(i, order, origin) for i = 0:order]

        tabfct = [LagrangePolynomial(convert.(T, coefficients)) for coefficients in tabfct_rat]
        return new{T,edge,order}(tabfct)

    end

end

function _c(k, n)
    coefficients = _getlagrangecoefficients(k, n, 0)
    result = zero(Rational{BigInt})
    for m in eachindex(coefficients)
        result += (isodd(m) ? coefficients[m] : -coefficients[m]) / m
    end
    return result
end
struct ABcoef
    tab::Array{Rational{Int}}
    function ABcoef(ordermax::Int)
        tab = zeros(Rational{Int}, ordermax, ordermax)
        for j = 1:ordermax, i = 1:j
            tab[i, j] = _c(i - 1, j - 1)
        end
        return new(tab)
    end
end
c(st::ABcoef, k, n) = st.tab[k, n]
