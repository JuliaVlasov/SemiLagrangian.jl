using Polynomials
import Base: +, *, -, ==, getindex, setindex!
# import Base: +, *
include("lapack.jl")

struct Spline{N}
    tabpol::Vector{Polynomials.Polynomial{Rational{N}}}
    function Spline(tabpol::Vector{Polynomials.Polynomial{Rational{N}}}) where{N<:Signed}
        return new{N}(tabpol)
    end
end

function Base.getindex(sp::Spline{N}, index) where{N<:Signed}
    i = index+1
    return if 1 <= i <= size(sp.tabpol, 1)
        sp.tabpol[i]
    else 
        zero(Polynomials.Polynomial{Rational{N}})
    end
end
function Base.setindex!(sp::Spline{N}, pol::Polynomials.Polynomial{Rational{N}}, index) where{N<:Signed}
    sp.tabpol[index-1] = pol
end
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
    i = Int64((floor(x)))+1 # it's different from ceil(x)
    if 1 <= i <= size(f.tabpol,1)
        return f.tabpol[i](x)
    else
        return zero(x)
    end
end

