"""
    InterpolationType{T, iscirc, order}

Abstract supertype for all interpolation type
"""
abstract type InterpolationType{T, iscirc, order} end
# modulo for "begin to one" array
modone(ind, n)=(n+ind-1)%n+1
gettabmod(lg)=modone.(1:3lg, lg)#
function get_allprecal(interp::InterpolationType{T,iscirc,order}, decint::Int, decfloat::T) where {T,iscirc,order}
    origin = -div(order, 2)
    inddeb , indfin = if decint < 0
        origin+decint, origin+order+decint
    else
        origin+decint, origin+order+decint
    end
    return [get_precal(interp, decfloat+i) for i=inddeb:indfin]
end
# get_allprecal(interp::InterpolationType{T}, decf::T) where {T}=get_allprecal(interp,0,decf)
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
    lgp = length(allprecal)
    # decintbeg = decint > 0 ? 0 : decint
    # decintend = decint < 0 ? 0 : decint
    # borne1=-decintbeg-origin-decintend
    # borne2=lg-decintend+origin-decintbeg-1
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



    # for i=1:lg
    #     indbeg = max(1,i+origin+decintbeg)
    #     indend = min(lg,indbeg+order)
    #     indbeg = indend-order
    #     ind = (i < lg/2) ? ind=i-indbeg+1 : lgp - (indend-i)
    # end 
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
function interpolate!( fp, fi, dec, interp::InterpolationType{T,iscirc, order}) where {T, iscirc, order}
    decint = convert(Int, floor(dec)) 
    decfloat = dec - decint
    if iscirc
        return interpolate!(fp, fi, decint, get_precal(interp, decfloat), interp )
    else
        return interpolate!(fp, fi, decint, get_allprecal(interp, decint, decfloat), interp )
    end
    missing
end


