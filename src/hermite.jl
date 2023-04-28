
using Polynomials

function _L(i, ord)
    ord % 2 == 1 || throw(DomainError("ord=$ord must be odd"))
    d = div(ord, 2)
    result = Polynomials.Polynomial([big(1 // 1)])
    for j = (-d):(d+1)
        if j != i
            result *= Polynomials.Polynomial([-j, 1 // 1] .// (i - j))
        end
    end
    return result
end

function _Lprim(i, ord)
    result = big(0 // 1)
    d = div(ord, 2)
    for j = (-d):(d+1)
        if i != j
            result += 1 // (i - j)
        end
    end
    return result
end
_K(i, ord) = _L(i, ord)^2 * Polynomial([-i, 1 // 1])
_H(i, ord) = _L(i, ord)^2 * (1 - 2 * _Lprim(i, ord) * Polynomial([-i, 1 // 1]))

function _bplus(i, rplus, splus)
    res = prod([big(-j // 1) for j = rplus:splus if j != 0 && j != i])
    return res * prod([big(1 // (i - j)) for j = rplus:splus if j != i])
end

function _bminus(i, rminus, sminus)
    return -_bplus(-i, -sminus, -rminus)
end


"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct PrecalHermite{ord,d}
    L::Vector{Polynomial{Rational{BigInt}}}
    Lprim::Vector{Rational{BigInt}}
    K::Vector{Polynomial{Rational{BigInt}}}
    H::Vector{Polynomial{Rational{BigInt}}}
    bplus::Vector{Rational{BigInt}}
    bminus::Vector{Rational{BigInt}}
    rplus::Int
    splus::Int
    rminus::Int
    sminus::Int
    function PrecalHermite(ord; flbis = false)
        ord % 2 == 1 || throw(DomainError("ord=$ord must be odd"))
        d = div(ord, 2)
        L = map(i -> _L(i, ord), (-d):(d+1))
        Lprim = map(i -> _Lprim(i, ord), (-d):(d+1))
        K = map(i -> _K(i, ord), (-d):(d+1))
        H = map(i -> _H(i, ord), (-d):(d+1))
        if flbis
            rplus = -d - 1
            splus = d
        else
            rplus = -d
            splus = d + 1
        end
        rminus = -splus
        sminus = -rplus

        # ref        bplus = map(i-> i != 0 ? _bplus(i,ord) : 0 , -d:d+1)
        bplus = map(i -> i != 0 ? _bplus(i, rplus, splus) : 0, rplus:splus)
        # ref        bminus = map(i-> i != 0 ? _bminus(i,ord) : 0 , -d-1:d)
        bminus = map(i -> i != 0 ? _bminus(i, rminus, sminus) : 0, rminus:sminus)
        bplus[1-rplus] = -sum(bplus)
        bminus[1-rminus] = -sum(bminus)
        return new{ord,d}(L, Lprim, K, H, bplus, bminus, rplus, splus, rminus, sminus)
    end
end

L(ph::PrecalHermite{ord,d}, i) where {ord,d} = ph.L[d+1+i]

Lprim(ph::PrecalHermite{ord,d}, i) where {ord,d} = ph.Lprim[d+1+i]

K(ph::PrecalHermite{ord,d}, i) where {ord,d} = ph.K[d+1+i]

H(ph::PrecalHermite{ord,d}, i) where {ord,d} = ph.H[d+1+i]

bplus(ph::PrecalHermite{ord,d}, i) where {ord,d} = ph.bplus[1-ph.rplus+i]

bminus(ph::PrecalHermite{ord,d}, i) where {ord,d} = ph.bminus[1-ph.rminus+i]

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct Hermite{T,edge,order} <: AbstractInterpolation{T,edge,order}
    tabfct::Vector{Polynomial{T}}
    function Hermite(
        order::Int,
        T::DataType = Float64;
        edge::EdgeType = CircEdge,
        flbis = false,
    )
        if flbis
            order % 4 == 3 || throw(DomainError("order=$order modulo 4 must equal to 3"))
            ord = div(order, 2)
        else
            order % 4 == 1 || throw(DomainError("order=$order modulo 4 must equal to 1"))
            ord = div(order, 2) + 1
        end
        decal = div(order + 1, 2)
        d = div(ord, 2)
        ph = PrecalHermite(ord; flbis = flbis)
        tabfct = map(x -> Polynomial(Rational{BigInt}[0]), 1:(order+1))
        for i = (-d):(d+1)
            tabfct[decal+i] += H(ph, i)
        end
        for i = (-d):0
            for k = (ph.rplus):(ph.splus)
                tabfct[decal+i+k] += bplus(ph, k) * K(ph, i)
            end
        end
        for i = 1:(d+1)
            for k = (ph.rminus):(ph.sminus)
                tabfct[decal+i+k] += bminus(ph, k) * K(ph, i)
            end
        end
        return new{T,edge,order}(tabfct)
    end
end
