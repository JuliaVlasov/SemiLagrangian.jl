
"""
$(TYPEDEF)

    BSpline{T, edge, order} <: AbstractInterpolation{T, edge, order}

Abstract supertype for all bspline interpolation type
"""
abstract type BSpline{T,edge,order} <: AbstractInterpolation{T,edge,order} end
isbspline(_::BSpline) = true
issolidentity(_::BSpline) = false
function get_kl_ku(order)
    ku = div(order, 2)
    kl = order - 1 - ku
    return kl, ku # kl <= ku
end
