
include("spline.jl")
abstract type InterpolationType{T, iscirc} end
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
    return res
end
@inline function get_precal!(v::Vector{T}, bsp::B_Spline{T}, decf::T) where{T}
    bs = get_bspline(bsp)
    @inbounds v .= [ bs[i](decf+i) for i=get_order(bsp):-1:0]
end
