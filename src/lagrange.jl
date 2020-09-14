using Polynomials
# using DynamicPolynomials

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
function _getpolylagrange(k::Int64, order::Int64, origin::Int64, N::DataType)
    0 <= k <= order || throw(DomainError("the constaint 0 <= k <= order is false"))
    N <: Union{BigInt,Int64} || throw(DomainError(N, "N must be Int64 or BigInt"))
    result = Polynomials.Polynomial([one(Rational{N})])
    for l=0:order
        if l != k
            result *= Polynomials.Polynomial([-(l+origin)//1,1//1])/(k-l)
        end
    end
    return result
end

function get_origin(order)
    v, _ = get_kl_ku(order)
    return -(v+1)
end

"""
    Lagrange{T, iscirc, granularity}
Lagrange Polynomials coefficients

# Fields :
- `coef::Matrix{T}` : Matrix with all Lagrange polynomials coefficients, each culumn from 0 to order had the coefficients for the corresponding Lagrange polynomial. the matrix had a size of (order + 1, order + 1).

"""
struct Lagrange{T, iscirc, granularity} <: InterpolationType{T, iscirc}
    coef::Matrix{T}
    function Lagrange(T::DataType, order; iscirc::Bool=true, granularity=1)
        type = order <= 10 ? Int64 : BigInt 
        coef = zeros(T, order+1, order+1)
        origin = get_origin(order)
        for i = 0:order
            coef[:,i+1] .= convert.(T, coeffs(_getpolylagrange( i, order, origin, type)))
        end
        new{T, iscirc, granularity}(coef) 
    end
    Lagrange(order; kwargs...)= Lagrange(Float64, order; kwargs...)
end
get_order(lag::Lagrange)= size(lag.coef,1)-1

"""
    polinterpol(
    lag::Lagrange, 
    resfct::Vector{T}
) where {T}
return the interpolation polynomial for the given values of a function

# Arguments
- `lag::Lagrange` : object with Lagrange coefficients
- `resfct::Vector{T}`` : result of functions for values lag.origin to lag.origin+size(lag.coef+1, 1)

# Returns :
- Polynomial{T} : the interpolation polynomial
"""
function polinterpol(
    lag::Lagrange, 
    resfct::Vector{T}
) where {T}
    return Polynomials.Polynomial(lag.coef*resfct)
end

# modulo for "begin to one" array
modone(ind, n)=(n+ind-1)%n+1
"""
    polinterpol(
    lag::Lagrange, 
    resfct::Vector{T},
    ind
) where {T<:Union{AbstractFloat,Complex{AbstractFloat}}}
return the interpolation polynomial for the given values of a function at a specified index

# Arguments
- `lag::Lagrange` : object with Lagrange coefficients
- `resfct::Vector{T}`` : result of functions for all values.
- `ind` : indices to take values from lag.origin +ind to ind + lag.origin+size(lag.coef+1, 1)

# Returns :
- Polynomial{T} : the interpolation polynomial
"""
function polinterpol(
    lag::Lagrange{T, iscirc, granularity}, 
    resfct::Vector{T},
    ind
) where {T, iscirc, granularity}
    order = get_order(lag)
    origin = get_origin(order)
    indbegin = origin+ind
    indend = indbegin+order
    listind = 
    if iscirc
        modone.(indbegin:indend, size(resfct,1))
    else
        decl = 0
        if indbegin < 1
            decl = indbegin - 1
            1:indend-decl
        elseif indend > size(resfct,1)
            decl = indend - size(resfct,1)
            indbegin-decl:size(resfct,1)
        else
            indbegin:indend
        end
    end
    f = resfct[listind]
    polret = polinterpol(lag, f)
    return if iscirc || decl == 0
        polret
    else
        polret(Polynomials.Polynomial([decl, one(T)]))
    end
end
"""
    interpolate!( fp, fi, dec, lag::Lagrange)
return the interpolation polynomial for the given values of a function a a specified index

# Arguments
- `fp` : output array of length n
- `fi` : input array of length n
- `dec` : offset in units of dx

# Returns :
- No return
"""
function interpolate!( adv, fp, fi, dec, 
    lag::Lagrange{T, iscirc, granularity}
) where {T, iscirc, granularity}
    if (dec >= 1 || dec < 1)
        cor = convert(Int, floor(dec))
        val = dec - cor
    else
        cor = 0
        val = dec
    end
    gr2, gr1=get_kl_ku(granularity)
    borne = size(fp,1)-gr2
    for i=gr1+1:granularity:borne
        pol = polinterpol(lag, fi, i+cor)
        for j=-gr1:gr2
            fp[i+j] = pol(j+val)
        end
    end
    reste = size(fp,1)%granularity
    if reste != 0
        gr1 = div(reste, 2)
        gr2 = reste-gr1-1
        i = size(fp,1)-gr2
        pol = polinterpol(lag, fi, i+cor)
        for j=-gr1:gr2
            fp[i+j] = pol(j+val)
        end
    end
end
get_order(lag::Lagrange{T}) where{T}=size(lag.coef,1)-1
get_type(lag::Lagrange{T, isc, gr}) where{T,isc,gr}="Lagrange{$T, $isc,$gr}"
