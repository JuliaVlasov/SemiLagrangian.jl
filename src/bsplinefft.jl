



struct B_SplineFFT{T, order, N} <: B_Spline{T, CircEdge, order}
    c_fft::Vector{Complex{T}}
    parfft
    fact_order::N
    tabpol::Vector{Polynomial{N}}
    function B_SplineFFT( order::Int, n::Int, T::DataType=Float64)
        bspline = SplineInt(order)
        tabpol = map(x -> bspline[order-x](Polynomial([order-x,1])), 0:order)
        N = typeof(bspline.fact_order)
        kl, ku = get_kl_ku(order)
        c=zeros(T,n)
        parfft=PrepareFftBig(n, T)
        tab_coef = convert.(T,bspline.(1:order))
        dec = n-kl-1
        for i=1:order
            c[(dec+i)%n+1] = tab_coef[i]
        end
        c_fft = fftgen(parfft, c)
        return new{T, order, N}(c_fft, parfft, bspline.fact_order, tabpol)
    end
    B_SplineFFT(o::Int, n::Int, _::T) where {T}=B_SplineFFT(o, n, T)
end
function sol(bsp::B_SplineFFT{T}, b::AbstractVector{T}) where{T}
    return real(ifftgen(bsp.parfft,fftgen(bsp.parfft,b) ./ bsp.c_fft))
end
get_n(bsp::B_SplineFFT) where{T}=size(bsp.c_fft,1)
get_order(bsp::B_SplineFFT{T, order}) where{T, order}=order
get_bspline(bsp::B_SplineFFT) where{T}=bsp.bspline
get_type(bsp::B_SplineFFT{T, order}) where{T, order}="B_SplineFFT{$T, $order}"
Base.show(io::IO, bsp::B_SplineFFT)=print(io, get_type(bsp))
istrace(bsp::B_SplineFFT)=false

get_fact_order(bsp::B_SplineFFT)=bsp.fact_order
get_tabpol(bsp::B_SplineFFT)=bsp.tabpol


