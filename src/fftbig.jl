#=
fft:
- Julia version: 
- Author: ymocquar
- Date: 2019-11-15
=#
# function that reverse the order of the pos lowest bits
using FFTW
function _reverse_num(num, pos)
    result = 0
    pos_m1 = pos-1
    for i=0:pos_m1
        if (num & (1 << i))  != 0
           result |= 1 << (pos_m1 - i)
        end
    end
    return result
end

"""
    PrepareFftBig( size_fft::Unsigned, [T=BigFloat])

Immutable structure to operate fft transform, 
x is the type of non transformed data also called signal.

# Arguments :
- `size_fft::Integer` : Number of values, must be a power of two
- `[T=BigFloat | x::T ]` : type of the values

# Implementation
- size_fft : size of the signal
- tab_permut : permutation
- root_one : size order roots of one
- root_one_conj : conjugate of root_one

"""

struct PrepareFftBig{T, NDIMS, NUMDIM}
    size_fft
    tab_permut
    root_one
    root_one_conj
    function PrepareFftBig( size_fft::Integer, x::T; ndims=2, numdim=1 ) where {T<:AbstractFloat}
        @assert prevpow(2,size_fft) == size_fft "size_fft=$size_fft is not a power of 2"
        power = convert(Int64,log2(size_fft))
        tab_permut = zeros(Int64,size_fft)
        root_one = zeros(Complex{T}, size_fft )
        root_one_conj = zeros(Complex{T}, size_fft )
        prec= precision(T)
        setprecision(prec+32) do
            for i=1:size_fft
                tab_permut[i] = _reverse_num( i-1, power) + 1
                root_one[i] = round(
    exp(one(T)*2big(pi)*im*(i-1)/size_fft),    
    digits=prec+16, 
    base=2 
)
                root_one_conj[i] = round(
    exp(-one(T)*2big(pi)*im*(i-1)/size_fft),
    digits=prec+16,
    base=2
)
            end
        end
        return new{T, ndims, numdim}(
    size_fft, 
    tab_permut, 
    Complex{T}.(root_one), 
    Complex{T}.(root_one_conj)
)
    end
end
function PrepareFftBig( size_fft::Integer, type::DataType; kwargs... ) 
    return PrepareFftBig( size_fft, one(type); kwargs...)
end
function PrepareFftBig( size_fft::Integer; kwargs...) 
    return PrepareFftBig( size_fft, one(BigFloat); kwargs...)
end

# Amazing that such function doesn't already exist
function permutselecteddim(in::Array{T}, numdim, perm) where {T}
    sz = size(in)
    out = zeros(T,sz)
    for i=1:sz[numdim]
        s_in = selectdim(in, numdim, i)
        s_out = selectdim(out, numdim, perm[i])
        s_out .= s_in
    end
    return out
end
function fftbig!(
    par::PrepareFftBig{T, NDIMS, NUMDIM}, 
    signal; 
    flag_inv=false
) where{T, NDIMS, NUMDIM}
    s=size(signal, NUMDIM)
    @assert prevpow(2,s) == s "size_fft(signal)=$s is not a power of 2"
    s_div2 = div(s,2)
    len = s
    n_len = len>>1
    nb_r = 1
    rootO = flag_inv ? par.root_one : par.root_one_conj;
#    prec= precision(real(rootO[1]))
#    setprecision(prec+32) do
         while n_len != 0
            start = 1
            suite = start+n_len
            for i=1:nb_r
                 deb = 1
                for j=deb:nb_r:s_div2
                    if NDIMS == 1
                        signal[start], signal[suite] = (signal[start] + signal[suite]), 
                        (signal[start] - signal[suite])*rootO[j]
                    else
                        s_start = selectdim(signal, NUMDIM, start)
                        s_suite = selectdim(signal, NUMDIM, suite)
                        s_sum = s_start+s_suite
                        s_diff = (s_start-s_suite)*rootO[j]
                        s_start .= s_sum
                        s_suite .= s_diff
                    end
                    # if NDIMS == 1
                    #     signal[start], signal[suite] = (signal[start] + signal[suite]), 
                    #     (signal[start] - signal[suite])*rootO[j]
                    # elseif NUMDIM == 1
                    #     signal[start, :], signal[suite, :] = (signal[start, :] + signal[suite, :]), 
                    #     (signal[start, :] - signal[suite, :])*rootO[j]
                    # else
                    #     signal[:,start], signal[:,suite] = (signal[:,start] + signal[:,suite]), 
                    #     (signal[:,start] - signal[:,suite])*rootO[j]
                    # end
                    start += 1
                    suite += 1
                end
                start = suite
                suite = start+n_len
            end
            len = n_len
            n_len >>= 1
            nb_r <<= 1
        end
#    end
# signal .= 
# if NDIMS == 1 
#     signal[par.tab_permut]
# elseif NUMDIM == 1
#     signal[par.tab_permut, :]
# else 
#     signal[:,par.tab_permut]
# end       
    signal .= 
    if NDIMS == 1 
        signal[par.tab_permut]
    else
        permutselecteddim(signal, NUMDIM, par.tab_permut)
    end       

    if flag_inv
        signal ./= s
    end
    return signal
end
function fftbig(
    par::PrepareFftBig{T, NDIMS, NUMDIM},
    signal;
    flag_inv=false
) where{T, NDIMS, NUMDIM}
    fl = flag_inv
    return fftbig!(
        par::PrepareFftBig,
        copy(convert(Array{Complex{T}}, signal)),
        flag_inv=fl
    )
end
fftgen(_::Any, t::Array{Complex{Float64}}) = fft(t, (1,))
fftgen!(_::Any, t::Array{Complex{Float64}}) = fft!(t, (1,))
fftgen(_::Any, t::Array{Float64}) = fft(t, (1,))
function fftgen(
    _::PrepareFftBig{T, NDIMS, NUMDIM}, 
    t::Array{Complex{Float64}}
) where {T, NDIMS, NUMDIM}
    return fft(t, (NUMDIM,))
end
function fftgen!(
    _::PrepareFftBig{T, NDIMS, NUMDIM}, 
    t::Array{Complex{Float64}}
) where {T, NDIMS, NUMDIM}
    return fft!(t, (NUMDIM,))
end
function fftgen(
    _::PrepareFftBig{T, NDIMS, NUMDIM}, 
    t::Array{Float64}
) where {T, NDIMS, NUMDIM}
    return fft(t, (NUMDIM,))
end
fftgen(p::PrepareFftBig, t::Array{T}) where {T<:AbstractFloat} = fftbig(p, t)
fftgen(p::PrepareFftBig, t::Array{Complex{T}}) where {T<:AbstractFloat} = fftbig(p, t)
fftgen!(p::PrepareFftBig, t::Array{Complex{T}}) where {T<:AbstractFloat} = fftbig!(p, t)
ifftgen(_::Any, t::Array{Complex{Float64}}) = ifft(t, (1,))
ifftgen!(_::Any, t::Array{Complex{Float64}}) = ifft!(t, (1,))
ifftgen(_::Any, t::Array{Float64}) = ifft(t, (1,))
function ifftgen(
    _::PrepareFftBig{T, NDIMS, NUMDIM}, 
    t::Array{Complex{Float64}}
) where {T, NDIMS, NUMDIM}
    return ifft(t, (NUMDIM,))
end
function ifftgen!(
    _::PrepareFftBig{T, NDIMS, NUMDIM}, 
    t::Array{Complex{Float64}}
) where {T, NDIMS, NUMDIM}
    return ifft!(t, (NUMDIM,))
end
function ifftgen(
    _::PrepareFftBig{T, NDIMS, NUMDIM}, 
    t::Array{Float64}
) where {T, NDIMS, NUMDIM}
    return ifft(t, (NUMDIM,))
end
ifftgen(p::PrepareFftBig, t::Array{T}) where {T}  = fftbig(p, t, flag_inv = true)
ifftgen(p::PrepareFftBig, t::Array{Complex{T}}) where {T}  = fftbig(p, t, flag_inv = true)
ifftgen!(p::PrepareFftBig, t::Array{Complex{T}}) where {T}  = fftbig!(p, t, flag_inv = true)
