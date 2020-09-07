
include("spline.jl")
abstract type B_Spline{T,iscirc} <: InterpolationType{T,iscirc} end
function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku
end

function interpolate!( adv, fp, fi, dec, 
    bsp::B_Spline{T, iscirc}
) where {T, iscirc}
    istrace(bsp) && println("begin interpolate with $(get_type(bsp))")
    res = sol(bsp, fi)
    n = get_n(bsp)
    order = get_order(bsp)
    kl, ku = get_kl_ku(order)

    decint = convert(Int, floor(dec))
    decfloat = dec-decint
    if decfloat > 0.5
        decfloat -= 1.
        decint += 1
    end

#    println("n=$n size(res)=$(size(res,1)) size(fi)=$(size(fi,1)) size(fp)=$(size(fp,1))")
    precal = get_bspline(bsp).((1:order) .- decfloat)
    if iscirc
        decall = n-kl+decint-2
        for i=1:n
            fp[i] = 0
            dec = decall + i
            for j=1:order
                # fp[i] += res[(i+n-bsp.ls.kl-2+j)%n+1]*bsp.bspline(j-1+dec)
                fp[i] += res[(dec+j)%n+1]*precal[j]
            end
            # diff = fp[i] -fi[i]
            # println("i2=$i diff=$diff")
        end
    else
        decall = -kl+decint-1
        for i=1:n
            # deb = max(1, bsp.ls.kl+2-i)
            # fin = min(order, n-i + bsp.ls.ku+1)
            fp[i] = 0
            dec = decall+i
            for j=1:order
                ind = dec+j
                if 1 <= ind <= n
                    fp[i] += res[ind]*precal[j]
                end
            end
        end
    end
    istrace(bsp) && println("end interpolate with $(get_type(bsp))")
end
