using Polynomials
using MultivariatePolynomials
using DynamicPolynomials
const MP = MultivariatePolynomials
const DP = DynamicPolynomials
const P = Polynomials
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
# if ! @isdefined XVar
    const XVar=PolyVar{true}("X")
    const YVar=PolyVar{true}("Y")
# end
function _getpolylagrange2d(k1::Int64, k2::Int64, order::Int64, origin::Int64, N::DataType)
    0 <= k1 <= order || throw(DomainError("the constaint 0 <= k1 <= order is false"))
    0 <= k2 <= order || throw(DomainError("the constaint 0 <= k2 <= order is false"))
    N <: Union{BigInt,Int64} || throw(DomainError(N, "N must be Int64 or BigInt"))
    result = DynamicPolynomials.Monomial([XVar, YVar], [0, 0])  # that is one
    result *= one(Rational{N})
    x1y0 = DynamicPolynomials.Monomial([XVar, YVar], [1, 0])  # that is X
    x0y1 = DynamicPolynomials.Monomial([XVar, YVar], [0, 1])  # that is Y
    for l=0:order
        if l != k1
            result *= (x1y0-(l+origin)//1)/(k1-l)
        end
        if l != k2
            result *= (x0y1-(l+origin)//1)/(k2-l)
        end
    end
    return result
end
function allcoefficients(
    p::DynamicPolynomials.Polynomial{true,T}, 
    ref::DynamicPolynomials.Polynomial{true}
) where {T}
    ind=1
    res = zeros(T, length(ref))
    a = DynamicPolynomials.coefficients(p)
    x = DynamicPolynomials.monomials(p)
    ref_x = DynamicPolynomials.monomials(ref)
    for i=1:length(ref)
#        println("i=$i ind=$ind")
        if x[ind] == ref_x[i]
            res[i] = a[ind]
            ind += 1
            if ind > length(p)
                break
            end
        end
    end
    return res
end
function polynomialfromall( a::Vector{T}, ref::DynamicPolynomials.Polynomial{true}) where {T}
    refmonovec = DynamicPolynomials.monomials(ref)
    z = collect(zip( a, refmonovec))
    z = filter( x-> x[1] != 0, z )
    return sum([ t[1]*t[2] for t in z])
end
"""
    Lagrange2d{iscirc, T, origin, granularity}
Lagrange Polynomials coefficients

# Fields :
- `coef::Matrix{Rational{N}}` : Matrice with all Lagrange polynomials coefficients, each culumn from 0 to order had the coefficients for the corresponding Lagrange polynomial. the matrix had a size of (order + 1, order + 1).
- `origin::Int64` : origin of the coefficients
- `iscirc::Bool` : 
"""
struct Lagrange2d{iscirc, T, origin, granularity} <: InterpolationType
    coef::Matrix{T}
    ref # polynomial with all possible monomials
    indref::Array{Int64,2} # correspond from x,y indices to ref index
    order
    function Lagrange2d(T::DataType, order; iscirc::Bool=true, granularity=1)
        origin = -div(order,2) 
        type = order <= 10 ? Int64 : BigInt
        op1 = order+1 
        len = op1^2
        coef = zeros(T, len, len)
        indref = zeros(Int64, op1, op1)
        ref = DynamicPolynomials.Monomial{true}()
        for i = 0:order, j = 0:order
            ref += DynamicPolynomials.Monomial([XVar, YVar], [i,j])
        end
        for (i, x) in enumerate(DynamicPolynomials.monomials(ref))
            expx = MP.exponents(x)
            indref[expx[1]+1,expx[2]+1] = i
        end
        for i = 0:order
            for j = 0:order
                indice = i*op1+j+1
                coef[:,indice] .= convert.(T, allcoefficients(_getpolylagrange2d( i, j, order, origin, type), ref))
            end
        end
        new{iscirc, T, origin, granularity}(coef, ref, indref, order) 
    end
    Lagrange2d(order; kwargs...)= Lagrange2d(Float64, order; kwargs...)
end
"""
    polinterpol(
    lag::Lagrange2d, 
    resfct::Vector{T}
) where {T<:Union{AbstractFloat,Complex{AbstractFloat}}}
return the interpolation polynomial for the given values of a function

# Arguments
- `lag::Lagrange` : object with Lagrange coefficients
- `resfct::Vector{T}`` : result of functions for values lag.origin to lag.origin+size(lag.coef+1, 1)

# Returns :
- Polynomial{T} : the interpolation polynomial
"""
function polinterpol(
    lag::Lagrange2d, 
    resfct::Vector{T}
) where {T<:Union{AbstractFloat,Complex{AbstractFloat}}}
    return polynomialfromall(lag.coef*resfct, lag.ref)
end

# modulo for "begin to one" array
modone(ind, n)=(n+ind-1)%n+1

function calindices(
    lag::Lagrange2d{iscirc, N, origin, granularity},
    ind::Int64,
    n::Int64
)where {iscirc, N, origin, granularity}
    decl = 0
    indbegin = origin+ind
    indend = indbegin+lag.order
    listind = 
    if iscirc
        modone.(indbegin:indend, n)
    else
        if indbegin < 1
            decl = indbegin - 1
            1:indend-decl
        elseif indend > n
            decl = indend - n
            indbegin-decl:n
        else
            indbegin:indend
        end
    end
    return listind, decl
end
"""
    polinterpol(
    lag::Lagrange2d, 
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
    lag::Lagrange2d{iscirc, N, origin, granularity}, 
    resfct::Array{T,2},
    ind::Tuple{Int64,Int64}
) where {iscirc, N, origin, T<:AbstractFloat, granularity}
    list_x, dec_x = calindices(lag, ind[1], size(resfct,1))
    list_y, dec_y = calindices(lag, ind[2], size(resfct,2))
    op1 = lag.order+1
    f = zeros(T,op1^2)
    for (i_x, x) in enumerate(list_x), (i_y, y) in enumerate(list_y)
        indice = (i_x-1)*op1+i_y
        f[indice] = resfct[x,y]
    end
    polret = polinterpol(lag, f)
    return if dec_x == 0 && dec_y == 0
        polret
    else
        x1y0 = DynamicPolynomials.Monomial([XVar, YVar], [1, 0])  # that is X
        x0y1 = DynamicPolynomials.Monomial([XVar, YVar], [0, 1])  # that is Y
        polret(x1y0+dec_x, x0y1+dec_y)
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
function interpolate!( adv, fp, fi, dec::Array{Tuple{T,T},2}, 
    lag::Lagrange2d{iscirc, N, origin, granularity}
) where {iscirc, N, origin, T<:AbstractFloat, granularity}
    gr1_x=div(granularity,2)
    gr2_x=granularity-gr1_x-1
    borne_x = size(fp,1)-gr2_x
    gr1_y=div(granularity,2)
    gr2_y=granularity-gr1_y-1
    borne_y = size(fp,2)-gr2_y
    for i_x=gr1_x+1:granularity:borne_x, i_y=gr1_y+1:granularity:borne_y
        pol = polinterpol(lag, fi, (i_x,i_y))
        for j_x=-gr1_x:gr2_x, j_y=-gr1_y:gr2_y
            x = i_x+j_x
            y = i_y+j_y
            # if x == y
            #     println("decx=$(dec[x,y][1]) decy=$(dec[x,y][2])")
            # end
            fp[x, y] = pol(dec[x,y][1], dec[x,y][2])
        end
    end
    @assert size(fp,1)%granularity == 0 "size(fp,1)%granularity != 0"
    @assert size(fp,2)%granularity == 0 "size(fp,2)%granularity != 0"
    # TODO cas ou reste != size%granularity pour chaque dimension
    # if reste != 0
    #     gr1 = div(reste, 2)
    #     gr2 = reste-gr1-1
    #     i = size(fp,1)-gr2
    #     pol = polinterpol(lag, fi, i+cor)
    #     for j=-gr1:gr2
    #         fp[i+j] = pol(j+val)
    #     end
    # end
end
