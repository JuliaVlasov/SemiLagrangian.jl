
"""
    B_Spline{T,iscirc, order} <: InterpolationType{T,iscirc, order}

Abstract supertype for all bspline interpolation type
"""
abstract type B_Spline{T,iscirc, order} <: InterpolationType{T,iscirc, order} end
isbspline(_::B_Spline)=true
function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku
end
function get_precal(bsp::B_Spline, decf)
    bs = get_bspline(bsp)
#    cor = get_order(bsp)%2 == 0
    res = [bs[i](decf+i)/bs.fact_order for i=get_order(bsp):-1:0]
    return res
end
@inline function get_precal!(v::Vector{T}, bsp::B_Spline{T}, decf::T) where{T}
    @inbounds v .= get_precal(bsp,decf)
end
