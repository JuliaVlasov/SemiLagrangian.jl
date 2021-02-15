@enum EdgeType CircEdge=1 InsideEdge=2

"""
    InterpolationType{T, edge, order}

Abstract supertype for all interpolation type
"""
abstract type InterpolationType{T, edge, order} end
@inline function get_precal(interp::InterpolationType{T}, decf::T) where{T}
    return @inbounds [T(fct(decf)) for fct in get_tabpol(interp)] ./ get_fact_order(interp)
end

@inline function get_precal!(v::Vector{T}, bsp::InterpolationType{T}, decf::T) where{T}
    @inbounds v .= get_precal(bsp,decf)
end

# modulo for "begin to one" array
modone(ind, n)=(n+ind-1)%n+1
gettabmod(lg)=modone.(1:3lg, lg)#
function get_allprecal(interp::InterpolationType{T, InsideEdge,order}, decint::Int, decfloat::T) where {T,order}
    origin = -div(order, 2)
    indbeg = origin+decint
    indend = indbeg+order
    return [get_precal(interp, decfloat+i) for i=indbeg:indend]
end

"""
    interpolate!( fp, fi, decint, precal, interp)

apply an offset to the function fi interpolate by interp struct, the result is in fp vector,
decint and precal are precompute with get_precal method.

# Arguments
- `fp::AbstractVector` : output vector
- `fi::AbstractVector` : input vector
- `decint` : offset in units of dx
- `precal::Vector` : vector of length order+1 precompute with get_precal(interp, dec) (dec is the offset)
- `interp::InterpolationType{T, CircEdge, order}` : interpolation implementation
- `tabmod=gettabmod(length(fi))` : precompute for "begin at one" modulo

# Returns :
- No return
"""
function interpolate!( 
    fp::AbstractVector{T}, 
    fi::AbstractVector{T},
    decint::Int, 
    precal::Vector{T}, 
    interp::InterpolationType{T, CircEdge, order},
    tabmod=gettabmod(length(fi))
) where {T, order}
    res = sol(interp,fi)
    origin = -div(order,2)
    lg = length(fi)
    @inbounds for i=1:lg
        indbeg=i+origin+decint+lg
        indend=indbeg+order
        fp[i] = sum(res[tabmod[indbeg:indend]] .* precal)
    end
    missing  
end
function interpolate!( 
    fp::AbstractVector{T}, fi::AbstractVector{T}, decint::Int, 
    allprecal::Vector{Vector{T}}, 
    interp::InterpolationType{T, InsideEdge, order},
    tabmod=gettabmod(length(fi))
) where {T, order}
    res = sol(interp,fi)
    origin = -div(order,2)
    lg = length(fi)
    lgp = length(allprecal)
    borne1=-decint-origin
    borne2=lg-decint+origin-1
    for i=1:borne1
        indbeg=1
        indend=order+1
        ind=i
#        @show 1, i, ind, indbeg, indend, decint, lgp
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end
    for i=borne1+1:borne2
        indbeg = i-borne1
        indend = indbeg+order
        ind = borne1+1
#        @show 2, i, ind, indbeg, indend, decint, lgp
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end
    for i=borne2+1:lg
        indbeg = lg-order
        indend = lg
        ind = lgp-(lg-i)
#        @show 3, i, ind, indbeg, indend, decint, lgp
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end
    missing  
end

"""
    interpolate!( fp, fi, dec, interp)

apply the offset dec to the function fi interpolate by interp struct, the result is in fp Vector

# Arguments
- `fp` : output vector of length n
- `fi` : input vector of length n
- `dec` : offset in units of dx
- `interp::InterpolationType` : interpolation implementation

# Returns :
- No return
"""
function interpolate!( fp, fi, dec, interp::InterpolationType{T, edge, order}) where {T, edge, order}
    decint = convert(Int, floor(dec)) 
    decfloat = dec - decint
    if edge == CircEdge
        return interpolate!(fp, fi, decint, get_precal(interp, decfloat), interp )
    else
        return interpolate!(fp, fi, decint, get_allprecal(interp, decint, decfloat), interp )
    end
    missing
end


