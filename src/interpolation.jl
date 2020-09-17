# modulo for "begin to one" array
modone(ind, n)=(n+ind-1)%n+1
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
    res = sol(interp,fi)
#    println("size res=$(size(res))")
    order=get_order(interp)
    # TODO : pas tres normal le +isbspline, a analyser
    #  mais pour l'intant ca fonctionne
    origin=get_origin(order) +isbspline(interp)
    if iscirc
        precal = get_precal(interp, decfloat)
    else
        allprecal = [get_precal(interp, decfloat+i) for i=origin:(origin+order)]
        # if isbspline(interp)
        #     allprecal=allprecal[end:-1:1]
        # end
        # println("size allprecal=$(size(allprecal))")
        println("allprecal=$allprecal")
        println("allprecal=$(convert.(Array{Float64,1},allprecal))")
    end
    n = size(fi,1)
    if iscirc
        # println("trace iscirc=true")
        for i=1:n
            indbeg=i+origin
            indend=indbeg+order
            fp[i] = sum(res[modone.(indbeg:indend, n)] .* precal)
        end
    else
        # println("trace iscirc=false")
        for i=1:n
            indbeg = max(1,i+origin)
#            indbeg = max(1,i+origin-isbspline(interp))
            indend = min(n,indbeg+order)
            indbeg = indend-order
            ind=i-indbeg+1
            fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
        end
    end
    
end
