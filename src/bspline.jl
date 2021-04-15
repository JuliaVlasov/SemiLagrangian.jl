
"""
    B_Spline{T, edge, order} <: AbstractInterpolation{T, edge, order}

Abstract supertype for all bspline interpolation type
"""
abstract type B_Spline{T,edge,order} <: AbstractInterpolation{T,edge,order} end
isbspline(_::B_Spline) = true
issolidentity(_::B_Spline) = false
function get_kl_ku(order)
    ku = div(order, 2)
    kl = order - 1 - ku
    return kl, ku # kl <= ku
end
