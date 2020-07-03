using Polynomials
abstract type InterpolationType end
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
    result = Polynomial([one(Rational{N})])
    for l=0:order
        if l != k
            result *= Polynomial([-(l+origin)//1,1//1])/(k-l)
        end
    end
    return result
end
"""
    LagrangeNew{iscirc, T, origin, granularity}
Lagrange Polynomials coefficients

# Fields :
- `coef::Matrix{Rational{N}}` : Matrice with all Lagrange polynomials coefficients, each culumn from 0 to order had the coefficients for the corresponding Lagrange polynomial. the matrix had a size of (order + 1, order + 1).
- `origin::Int64` : origin of the coefficients
- `iscirc::Bool` : 
"""
struct LagrangeNew{iscirc, T, origin, granularity} <: InterpolationType
    coef::Matrix{T}
    function LagrangeNew(T::DataType, order; iscirc::Bool=true, granularity=1)
        origin = -div(order,2) 
        type = order <= 20 ? Int64 : BigInt 
        coef = zeros(T, order+1, order+1)
        for i = 0:order
            coef[:,i+1] .= convert.(T, coeffs(_getpolylagrange( i, order, origin, type)))
        end
        new{iscirc, T, origin, granularity}(coef) 
    end
    LagrangeNew(order; kwargs...)= LagrangeNew(Float64, order; kwargs...)
end
"""
    polinterpol(
    lag::LagrangeNew, 
    resfct::Vector{T}
) where {T<:Union{AbstractFloat,Complex{AbstractFloat}}}
return the interpolation polynomial for the given values of a function

# Arguments
- `lag::LagrangeNew` : object with Lagrange coefficients
- `resfct::Vector{T}`` : result of functions for values lag.origin to lag.origin+size(lag.coef+1, 1)

# Returns :
- Polynomial{T} : the interpolation polynomial
"""
function polinterpol(
    lag::LagrangeNew, 
    resfct::Vector{T}
) where {T<:Union{AbstractFloat,Complex{AbstractFloat}}}
    return Polynomial(lag.coef*resfct)
end

# modulo for "begin to one" array
modone(ind, n)=(n+ind-1)%n+1
"""
    polinterpol(
    lag::LagrangeNew, 
    resfct::Vector{T},
    ind
) where {T<:Union{AbstractFloat,Complex{AbstractFloat}}}
return the interpolation polynomial for the given values of a function at a specified index

# Arguments
- `lag::LagrangeNew` : object with Lagrange coefficients
- `resfct::Vector{T}`` : result of functions for all values.
- `ind` : indices to take values from lag.origin +ind to ind + lag.origin+size(lag.coef+1, 1)

# Returns :
- Polynomial{T} : the interpolation polynomial
"""
function polinterpol(
    lag::LagrangeNew{iscirc, N, origin, granularity}, 
    resfct::Vector{T},
    ind
) where {iscirc, N, origin, T<:AbstractFloat, granularity}
    indbegin = origin+ind
    indend = indbegin+size(lag.coef,1)-1
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
        polret(Polynomial([decl, one(T)]))
    end
end
"""
    interpolate!( fp, fi, dec, lag::LagrangeNew)
return the interpolation polynomial for the given values of a function a a specified index

# Arguments
- `fp` : output array of length n
- `fi` : input array of length n
- `dec` : offset in units of dx

# Returns :
- No return
"""
function interpolate!( adv, fp, fi, dec, 
    lag::LagrangeNew{iscirc, N, origin, granularity}
) where {iscirc, N, origin, T<:AbstractFloat, granularity}
    # if (dec >= 1 || dec < 1)
    #     cor = Int64(floor(dec))
    #     val = dec - cor
    # else
        cor = 0
        val = dec
    # end
    gr1=div(granularity,2)
    gr2=granularity-gr1-1
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
