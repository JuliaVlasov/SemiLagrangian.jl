
include("spline.jl")
abstract type B_Spline{T,iscirc} <: InterpolationType{T,iscirc} end
isbspline(_::B_Spline)=true
function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku
end
function get_precal(bsp, decf)
    bs = get_bspline(bsp)
    res = [ bs[i](decf+i) for i=get_order(bsp):-1:0]
#    res = [ bs[i](1+i - decf) for i=0:get_order(bsp)]
    @show res
    return res
end
# function interpolate!( adv, fp, fi, dec, 
#     bsp::B_Spline{T, iscirc}
# ) where {T, iscirc}
#     istrace(bsp) && println("begin interpolate with $(get_type(bsp))")
#     res = sol(bsp, fi)
#     # println("res=$res")
#     # println("res=$(convert.(Float64, res))")
#     n = get_n(bsp)
#     order = get_order(bsp)
#     kl, ku = get_kl_ku(order)

#     decint = convert(Int, floor(dec))
#     decfloat = dec-decint

#     println("decint=$decint decfloat=$decfloat")
#     # we always have 0 <= decfloat < 1

    
 
# #    println("n=$n size(res)=$(size(res,1)) size(fi)=$(size(fi,1)) size(fp)=$(size(fp,1))")
#     precal = get_bspline(bsp).(((order):-1:0) .+ decfloat)
# #    println("precal=$(convert.(Float64,precal)), $precal")
#     if iscirc
# #        println("trace iscirc=true")
#         decall = n-kl+decint-2
#         for i=1:n
#             fp[i] = 0
#             decloc = decall + i
#             for j=1:(order+1)
#                 # fp[i] += res[(i+n-bsp.ls.kl-2+j)%n+1]*bsp.bspline(j-1+dec)
#                 fp[i] += res[(decloc+j)%n+1]*precal[j]
#             end
#             # diff = fp[i] -fi[i]
#             # println("i2=$i diff=$diff")
#         end
#     else
# #        println("trace iscirc=false")
#         decall = -kl+decint-2
#         for i=1:n
#             # deb = max(1, bsp.ls.kl+2-i)
#             # fin = min(order, n-i + bsp.ls.ku+1)
#             fp[i] = 0
#             decloc = decall+i
#             for j=1:(order+1)
#                 ind = decloc+j
#                 if 1 <= ind <= n
#                     # if i <= 2 || i >= n-1
#                     #     println("ind=$ind j=$j")
#                     # end
#                     fp[i] += res[ind]*precal[j]
#                 end
#             end
#         end
#     end
#     istrace(bsp) && println("end interpolate with $(get_type(bsp))")
# end
