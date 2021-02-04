abstract type InterpolationType{T, iscirc, order} end
# modulo for "begin to one" array
modone(ind, n)=(n+ind-1)%n+1
gettabmod(lg)=modone.(1:3lg, lg)
function get_allprecal(interp::InterpolationType{T,iscirc,order}, decfloat::T) where {T,iscirc,order}
    origin = -div(order, 2)
    return [get_precal(interp, decfloat+i) for i=origin:(origin+order)]
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
- `interp::InterpolationType` : interpolation implementation

# Returns :
- No return
"""
function interpolate!( 
    fp::AbstractVector{T}, fi::AbstractVector{T}, decint::Int, 
    precal::Vector{T}, 
    interp::InterpolationType{T,true, order},
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
    interp::InterpolationType{T, false, order},
    tabmod=gettabmod(length(fi))
) where {T, order}
    res = sol(interp,fi)
    origin = -div(order,2)
    lg = length(fi)
    @inbounds for i=1:lg
        indbeg = max(1,i+origin)
        indend = min(lg,indbeg+order)
        indbeg = indend-order
        ind=i-indbeg+1
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end 
    missing  
end


function interpolate!( fp, fi, decint, precal::Vector{Vector{T}}, interp::InterpolationType{T,false, order}, tabmod=gettabmod(length(fi))) where {T, order}
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
function interpolate!( fp, fi, dec, interp::InterpolationType{T,iscirc}) where {T, iscirc}
    decint = convert(Int, floor(dec))
    decfloat = dec - decint
    precal = if iscirc
        get_precal(interp, decfloat)
    else
        get_allprecal(interp, decfloat)
    end
    if iscirc
        return interpolate!(fp, fi, decint, get_precal(interp, decfloat), interp )
    else
        return interpolate!(fp, fi, decint, get_allprecal(interp, decfloat), interp )
    end
    missing
end


