
"""
    B_Spline{T, edge, order} <: AbstractInterpolation{T, edge, order, 1}

Abstract supertype for all bspline interpolation type
"""
abstract type B_Spline{T, edge, order} <: AbstractInterpolation{T, edge, order, 1} end
isbspline(_::B_Spline)=true
function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku # kl <= ku
end

