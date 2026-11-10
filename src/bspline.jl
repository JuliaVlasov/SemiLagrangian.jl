"""
    BSpline{T, bc, order} <: AbstractInterpolation{T, bc, order}

Abstract supertype for all bspline interpolation type
"""
abstract type BSpline{T, bc, order} <: AbstractInterpolation{T, bc, order} end
isbspline(_::BSpline) = true
is_preprocess_identity(_::BSpline) = false
function get_kl_ku(order)
    ku = div(order, 2)
    kl = order - 1 - ku
    return kl, ku # kl <= ku
end
