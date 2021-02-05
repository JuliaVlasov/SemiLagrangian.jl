

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

# @inline function get_origin(order)
#     return -div(order,2)
# end

"""
    Lagrange{T, iscirc, order, N}
    Lagrange(T::DataType, order; iscirc::Bool=true)
Lagrange Polynomials coefficients

# Type parameters
- 'T' : the type of data that is interpolate
- 'iscirc::Bool' : true if function is circular
- `order::Int`: order of lagrange interpolation
- `N` : type of integer, in fact Int64 or BigInt that is used to store lagrange polynomial

# Implementation :
- `fact_order::N` : factorial of the order
- `lagpol:Vector{Polynomial{N}}` : vector of all lagrange polynomial, per example the k-th Lagrange polynomial for the order is lagpol[k+1]/fact_order

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
@inline get_precal(lag::Lagrange{T},decf) where{T}=@inbounds [T(fct(decf))/lag.fact_order for fct in lag.lagpol]
@inline get_precal!(v::Vector{T}, lag::Lagrange{T},decf) where{T}=@inbounds v .= get_precal(lag, decf)
@inline sol(lag::Lagrange,b)=b
@inline isbspline(_::Lagrange)=false
# function print(pol::Lagrange)
#     println("type=$(get_type(pol))")
#     println("polynomes=$(pol.lagpol)")
# end
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
