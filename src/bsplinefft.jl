

include("../src/fftbig.jl")

struct B_SplineFFT{T} <: B_Spline{T,true}
    order
    bspline::Spline
    c_fft::Vector{Complex{T}}
    parfft
    function B_SplineFFT( order, n, eltfortype::T) where{T}
        bspline = getbspline(order, 0)
        kl, ku = get_kl_ku(order)
        c=zeros(T,n)
        parfft=PrepareFftBig(n, T)
        tab_coef = convert.(T,bspline.(1:order))
        dec = n-kl-1
        for i=1:order
            c[(dec+i)%n+1] = tab_coef[i]
        end
        c_fft = fftgen(parfft, c)
        return new{T}(order, bspline, c_fft, parfft)
    end
    B_SplineFFT(o, n, t::DataType ; kwargs... )=B_SplineFFT(o, n, one(t) ; kwargs... )

end
function sol(bsp::B_SplineFFT{T}, b::Vector{T}) where{T}
    return real(ifftgen(bsp.parfft,fftgen(bsp.parfft,b) ./ bsp.c_fft))
end
get_n(bsp::B_SplineFFT) where{T}=size(bsp.c_fft,1)
get_order(bsp::B_SplineFFT) where{T}=bsp.order
get_bspline(bsp::B_SplineFFT) where{T}=bsp.bspline
get_type(bsp::B_SplineFFT{T}) where{T}="B_SplineFFT{$T}"
istrace(bsp::B_SplineFFT)=false

