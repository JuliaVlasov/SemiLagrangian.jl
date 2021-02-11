
"""
    B_Spline{T,iscirc, order} <: InterpolationType{T,iscirc, order}

Abstract supertype for all bspline interpolation type
"""
abstract type B_Spline{T,iscirc, order} <: InterpolationType{T,iscirc, order} end
isbspline(_::B_Spline)=true
function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku # kl <= ku
end
# function get_kl_ku(order)
#     kl = div(order,2)
#     ku = order-1-kl
#     return kl, ku # kl >= ku
# end
# function get_precal(bsp::B_Spline, decf)
#     bs = get_bspline(bsp)
# #    cor = get_order(bsp)%2 == 0
#     res = [bs[i](decf+i)/bs.fact_order for i=get_order(bsp):-1:0]
# #    res = [bs[0](i+1-decf)/bs.fact_order for i=0:get_order(bsp)]
#     return res
# end
@inline get_precal(bsp::B_Spline{T}, decf) where{T}= @inbounds [T(fct(decf)) for fct in get_tabpol(bsp)] ./ get_fact_order(bsp)

@inline function get_precal!(v::Vector{T}, bsp::B_Spline{T}, decf::T) where{T}
    @inbounds v .= get_precal(bsp,decf)
end
