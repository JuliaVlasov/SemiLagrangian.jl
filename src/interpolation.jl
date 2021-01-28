# modulo for "begin to one" array
modone(ind, n)=(n+ind-1)%n+1
gettabmod(lg)=modone.(1:3lg, lg)

function interpolate!( fp, fi, decint, precal, interp::InterpolationType{T,true}, tabmod=gettabmod(length(fi))) where {T}
    res = sol(interp,fi)
    order = size(precal,1)-1
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
return the interpolation polynomial for the given values of a function a a specified index

# Arguments
- `fp` : output array of length n
- `fi` : input array of length n
- `dec` : offset in units of dx

# Returns :
- No return
"""
function interpolate!( fp, fi, dec, interp::InterpolationType{T,iscirc}) where {T, iscirc}
    decint = convert(Int, floor(dec))
    decfloat = dec - decint
    if iscirc
        return interpolate!(fp, fi, decint, get_precal(interp, decfloat), interp )
    else
        res = sol(interp,fi)
        order = get_order(interp)
        origin = -div(order, 2)
        allprecal = [get_precal(interp, decfloat+i) for i=origin:(origin+order)]
        n = length(fi)
        for i=1:n
            indbeg = max(1,i+origin)
            indend = min(n,indbeg+order)
            indbeg = indend-order
            ind=i-indbeg+1
            fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
        end
    end
    missing
end


