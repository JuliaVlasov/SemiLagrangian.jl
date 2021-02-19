
import Base: +, *, -, ==, getindex, setindex!
# import Base: +, *

abstract type AbstractSpline{N} end
struct Spline{N} <: AbstractSpline{N}
    tabpol::Vector{Polynomials.Polynomial{Rational{N}}}
    function Spline(tabpol::Vector{Polynomials.Polynomial{Rational{N}}}) where{N<:Signed}
        return new{N}(tabpol)
    end
end

function Base.getindex(sp::AbstractSpline{N}, index::Integer) where{N<:Signed}
    i = index + 1
    return if 1 <= i <= size(sp.tabpol, 1)
        sp.tabpol[i]
    else 
        zero(Polynomials.Polynomial{Rational{N}})
    end
end
# function Base.getindex(sp::AbstractSpline{N}, index::AbstractRange) where{N<:Signed}
#     return Base.getindex.((sp::AbstractSpline{N},), index::AbstractRange)
# end
# function Base.setindex!(sp::AbstractSpline{N}, pol::Polynomials.Polynomial{Rational{N}}, index) where{N<:Signed}
#     sp.tabpol[index-1] = pol
# end
Base.size(sp::AbstractSpline, dim=1)=size(sp.tabpol,1)
function +(a::Spline{N}, b::Spline{N}) where{N<:Signed}
    sizenew = max(size(a.tabpol,1), size(b.tabpol,1))
    tabpolnew = zeros(Polynomials.Polynomial{Rational{N}},sizenew)
    for i=1:sizenew
        tabpolnew[i] += a[i-1]+b[i-1]
    end
    return(Spline(tabpolnew))
end
function -(a::Spline{N}, b::Spline{N}) where{N<:Signed}
    sizenew = max(size(a.tabpol,1), size(b.tabpol,1))
    tabpolnew = zeros(Polynomials.Polynomial{Rational{N}},sizenew)
    for i=1:sizenew
        tabpolnew[i] += a[i-1]-b[i-1]
    end
    return(Spline(tabpolnew))
end
function ==(a::Spline{N}, b::Spline{N}) where{N<:Signed}
    return a.tabpol == b.tabpol
end

function *( a::Spline{N}, pol::Polynomials.Polynomial{Rational{N}}) where{N<:Signed}
    tabpol2 = deepcopy(a.tabpol)
    for i=1:size(tabpol2,1)
        tabpol2[i] *= pol
    end
    return Spline(tabpol2)
end

function decal( a::Spline{N}, n) where{N<:Signed}
    return if n == 0
        a
    else
        tabpol = deepcopy(a.tabpol)
        poldec = Polynomials.Polynomial([-n//1,one(Rational{N})])
        for i=1:size(a.tabpol,1)
            tabpol[i] = a.tabpol[i](poldec)
        end
        Spline(vcat(zeros(Polynomials.Polynomial{Rational{N}},n),tabpol))
    end
end
w(p, j)=Polynomials.Polynomial([-j//p, 1//p])
function _getbspline(n::N, j)  where{N<:Signed}
 #   println("_getbspline($n=$n , j=$j ) N=$N")
    if n == zero(N) 
        ret = decal(Spline(ones(Polynomials.Polynomial{Rational{N}},1)), j)
    else
        n1 = _getbspline( n-1, j)
        n2 = decal(n1,1)
        ret = n1*w(n,j) +n2*(1-w(n,j+1))
    end
#    println("n=$n j=$j N=$N ret=$ret")
    return ret
end
function getbspline(n, j)
    return if ( n > 8 )
        _getbspline( big(n), j)
    else
        _getbspline( n, j)
    end
end
function (f::Spline{N})(x) where{N<:Signed}
    i = Int64((floor(x))) # it's different from ceil(x)
    if 0 <= i < size(f,1)
        return f[i](x)
    else
        return zero(x)
    end
end
struct SplineInt{N} <: AbstractSpline{N}
    fact_order::N
    tabpol::Vector{Polynomials.Polynomial{N}}
end
# function SplineInt(order)
#     sp = getbspline(order, 0)
#     N = order <= 13 ? Int64 : BigInt
#     fact_order = factorial(big(order))
#     return SplineInt{N}(N(fact_order), map( x->Polynomial(N.(fact_order*coeffs(sp[x]))), 0:order))
# end
# function (f::SplineInt{N})(x::T) where{N<:Signed, T <: AbstractFloat}
#     i = Int64((floor(x))) # it's different from ceil(x)
#     if 0 <= i < size(f,1)
#         return f[i](x)/f.fact_order
#     else
#         return zero(x)/f.fact_order
#     end
# end
# function (f::SplineInt{N})(x::Union{Rational,Int}) where{N<:Signed}
#     i = Int64((floor(x))) # it's different from ceil(x)
#     if 0 <= i < size(f,1)
#         return f[i](x)//f.fact_order
#     else
#         return zero(x)//f.fact_order
#     end
# end
