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

@inline function get_origin(order)
    v, _ = get_kl_ku(order)
    return -(v+1)
end

"""
    Lagrange{T, iscirc, granularity}
Lagrange Polynomials coefficients

# Fields :
- `coef::Matrix{T}` : Matrix with all Lagrange polynomials coefficients, each culumn from 0 to order had the coefficients for the corresponding Lagrange polynomial. the matrix had a size of (order + 1, order + 1).

"""
# struct Lagrange{T, iscirc, granularity} <: InterpolationType{T, iscirc}
#     coef::Matrix{T}
#     function Lagrange(T::DataType, order; iscirc::Bool=true, granularity=1)
#         type = order <= 10 ? Int64 : BigInt 
#         coef = zeros(T, order+1, order+1)
#         origin = get_origin(order)
#         for i = 0:order
#             coef[:,i+1] .= convert.(T, coeffs(_getpolylagrange( i, order, origin, type)))
#         end
#         new{T, iscirc, granularity}(coef) 
#     end
#     Lagrange(order; kwargs...)= Lagrange(Float64, order; kwargs...)
# end
# get_order(lag::Lagrange)= size(lag.coef,1)-1

struct Lagrange{T, iscirc} <: InterpolationType{T, iscirc}
    lagpol
    function Lagrange(T::DataType, order; iscirc::Bool=true) 
        type = order <= 10 ? Int64 : BigInt 
        origin = get_origin(order)
        lagpol = [_getpolylagrange( i, order, origin, type) for i=0:order]
        new{T, iscirc}(lagpol) 
    end
    Lagrange(order; kwargs...)= Lagrange(Float64, order; kwargs...)
end
@inline get_order(lag::Lagrange)= size(lag.lagpol,1)-1
@inline get_type(lag::Lagrange{T, isc}) where{T,isc}="Lagrange{$T, $isc}"
@inline get_precal(lag::Lagrange,decf)=[fct(decf) for fct in lag.lagpol]
@inline sol(lag::Lagrange,b)=b
@inline isbspline(_::Lagrange)=false
# """
#     polinterpol(
#     lag::Lagrange, 
#     resfct::Vector{T}
# ) where {T}
# return the interpolation polynomial for the given values of a function

# # Arguments
# - `lag::Lagrange` : object with Lagrange coefficients
# - `resfct::Vector{T}`` : result of functions for values lag.origin to lag.origin+size(lag.coef+1, 1)

# # Returns :
# - Polynomial{T} : the interpolation polynomial
# """
# function polinterpol(
#     lag::Lagrange, 
#     resfct::Vector{T}
# ) where {T}
#     return Polynomials.Polynomial(lag.coef*resfct)
# end

# # modulo for "begin to one" array
# modone(ind, n)=(n+ind-1)%n+1
# """
#     polinterpol(
#     lag::Lagrange, 
#     resfct::Vector{T},
#     ind
# ) where {T<:Union{AbstractFloat,Complex{AbstractFloat}}}
# return the interpolation polynomial for the given values of a function at a specified index

# # Arguments
# - `lag::Lagrange` : object with Lagrange coefficients
# - `resfct::Vector{T}`` : result of functions for all values.
# - `ind` : indices to take values from lag.origin +ind to ind + lag.origin+size(lag.coef+1, 1)

# # Returns :
# - Polynomial{T} : the interpolation polynomial
# """
# function polinterpol(
#     lag::Lagrange{T, iscirc, granularity}, 
#     resfct::Vector{T},
#     ind
# ) where {T, iscirc, granularity}
#     order = get_order(lag)
#     origin = get_origin(order)
#     indbegin = origin+ind
#     indend = indbegin+order
#     listind = 
#     if iscirc
#         modone.(indbegin:indend, size(resfct,1))
#     else
#         decl = 0
#         if indbegin < 1
#             decl = indbegin - 1
#             1:indend-decl
#         elseif indend > size(resfct,1)
#             decl = indend - size(resfct,1)
#             indbegin-decl:size(resfct,1)
#         else
#             indbegin:indend
#         end
#     end
#     f = resfct[listind]
#     polret = polinterpol(lag, f)
#     return if iscirc || decl == 0
#         polret
#     else
#         polret(Polynomials.Polynomial([decl, one(T)]))
#     end
# end
# modulo for "begin to one" array
# modone(ind, n)=(n+ind-1)%n+1
# """
#     interpolate!( fp, fi, dec, lag::Lagrange)
# return the interpolation polynomial for the given values of a function a a specified index

# # Arguments
# - `fp` : output array of length n
# - `fi` : input array of length n
# - `dec` : offset in units of dx

# # Returns :
# - No return
# """
# function interpolate!( adv, fp, fi, dec, 
#     lag::Lagrange{T, iscirc}
# ) where {T, iscirc, granularity}
#     decint = convert(Int, floor(dec))
#     decfloat = dec - decint

#     res = sol(lag,fi)

#     println("size res=$(size(res))")

#     order=get_order(lag)
#     origin=get_origin(order)

#     if iscirc
#         precal = get_precal(lag, decfloat)
#     else
#         allprecal = [get_precal(lag, decfloat+i) for i=origin:(origin+order)]
#         println("size allprecal=$(size(allprecal))")
#         println("allprecal=$allprecal")
#         println("allprecal=$(convert.(Array{Float64,1},allprecal))")
#     end


#     n = size(fi,1)
#     if iscirc
#         println("trace iscirc=true")
#         for i=1:n
#             indbeg=i+origin
#             indend=indbeg+order
#             fp[i] = res[modone.(indbeg:indend, n)] .* precal
#         end
#     else
#         println("trace iscirc=false")
#         for i=1:n
#             indbeg = max(1,i+origin)
#             indend = min(n,indbeg+order)
#             indbeg = indend-order
#             ind=i-indbeg+1
#             fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
#         end
#     end
    
# end
