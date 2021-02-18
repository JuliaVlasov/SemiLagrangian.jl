



struct B_SplineFFT{T, order} <: B_Spline{T, CircEdge, order}
    c_fft::Vector{Complex{T}}
    parfft::PrepareFftBig
    tabfct::Vector{Polynomial{T}}
    function B_SplineFFT( order::Int, n::Int, T::DataType=Float64)
        # bspline = SplineInt(order)
        # tabfct = map(x -> bspline[order-x](Polynomial([order-x,1])), 0:order)
        bspline = getbspline(order, 0)
        tabfct_rat = map(x -> bspline[order-x](Polynomial([order-x,1])), 0:order)
        kl, ku = get_kl_ku(order)
        c=zeros(T,n)
        parfft=PrepareFftBig(n, T)
        tab_coef = convert.(T,bspline.(1:order))
        dec = n-kl-1
        for i=1:order
            c[(dec+i)%n+1] = tab_coef[i]
        end
        c_fft = fftgen(parfft, c)
        return new{T, order}(c_fft, parfft, convert.(Polynomial{T}, tabfct_rat))
    end
    B_SplineFFT(o::Int, n::Int, _::T) where {T}=B_SplineFFT(o, n, T)
end
function sol(bsp::B_SplineFFT{T}, b::AbstractVector{T}) where{T}
    return real(ifftgen(bsp.parfft,fftgen(bsp.parfft,b) ./ bsp.c_fft))
end
get_n(bsp::B_SplineFFT) where{T}=size(bsp.c_fft,1)
# get_bspline(bsp::B_SplineFFT) where{T}=bsp.bspline


