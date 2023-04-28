
"""
$(TYPEDEF)

    BSplineFFT{T, order} <: AbstractInterpolation{T, CircEdge, order}
    BSplineFFT( order::Int, n::Int, T::DataType=Float64)

Type containing spline coefficients for b-spline interpolation based on fft, using the fact that b-spline matrix is a circulant matrix

# Type parameters
- `T` : the type of data that is interpolate
- `order::Int`: order of lagrange interpolation

# Implementation :
- `c_fft::Vector{Complex{T}}` : fft transform of coefficients
- `parfft::PrepareFftBig` : fft precomputed data
- `tabfct::Vector{Polynomial{T}}` : function table for interpolation

# Arguments : 
- `n` : size of the matrix
- `order` : the order of interpolation
- `[T::DataType=Float64]` : The type values to interpolate 

"""
struct BSplineFFT{T,order} <: BSpline{T,CircEdge,order}
    c_fft::Vector{Complex{T}}
    parfft::PrepareFftBig
    tabfct::Vector{Polynomial{T}}
    function BSplineFFT(order::Int, n::Int, T::DataType = Float64)
        bspline = getbspline(order, 0)
        tabfct_rat = map(x -> bspline[order-x](Polynomial([order - x, 1])), 0:order)
        kl, ku = get_kl_ku(order)
        c = zeros(T, n)
        parfft = PrepareFftBig(n, T)
        tab_coef = convert.(T, bspline.(1:order))
        dec = n - kl - 1
        for i = 1:order
            c[(dec+i)%n+1] = tab_coef[i]
        end
        c_fft = fftgen(parfft, c)
        return new{T,order}(c_fft, parfft, convert.(Polynomial{T}, tabfct_rat))
    end
    BSplineFFT(o::Int, n::Int, _::T) where {T} = BSplineFFT(o, n, T)
end

"""
$(SIGNATURES)
"""
function sol(bsp::BSplineFFT{T}, b::AbstractVector{T}) where {T}
    return real(ifftgen(bsp.parfft, fftgen(bsp.parfft, b) ./ bsp.c_fft))
end

"""
$(SIGNATURES)
"""
function sol!(Y::AbstractVector{T}, bsp::BSplineFFT{T}, b::AbstractVector{T}) where {T}
    return Y .= real(ifftgen(bsp.parfft, fftgen(bsp.parfft, b) ./ bsp.c_fft))
end
