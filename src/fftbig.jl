"""
$(SIGNATURES)

function that reverse the order of the pos lowest bits
"""
function _reverse_num(num, pos)
    result = 0
    pos_m1 = pos - 1
    for i = 0:pos_m1
        if (num & (1 << i)) != 0
            result |= 1 << (pos_m1 - i)
        end
    end
    return result
end

"""
$(TYPEDEF)

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

# struct parameters
- T : type of Float 
- NUMDIMS : Numbers of dimensions
- DIMS : tuple of dimensions

"""
struct PrepareFftBig{T,NUMDIMS,DIMS}
    size_fft::Any
    tab_permut::Any
    root_one::Any
    root_one_conj::Any
    function PrepareFftBig(
        size_fft::NTuple{N,Integer},
        x::T;
        numdims = 2,
        dims::NTuple{N,Integer} = (1,),
    ) where {T<:AbstractFloat,N}
        tab_permut = Vector{Vector{Int64}}(undef, N)
        root_one = Vector{Vector{Complex{T}}}(undef, N)
        root_one_conj = Vector{Vector{Complex{T}}}(undef, N)
        for (inddim, dim) in enumerate(dims)
            sfft = size_fft[inddim]
            @assert prevpow(2, sfft) == sfft "size_fft[$inddim]=$sfft is not a power of 2"
            power = convert(Int64, log2(sfft))
            prec = precision(T)
            setprecision(prec + 32) do
                tab_permut[inddim] = tperm = Vector{Int64}(undef, sfft)
                root_one[inddim] = rone = Vector{Complex{T}}(undef, sfft)
                root_one_conj[inddim] = rone_conj = Vector{Complex{T}}(undef, sfft)
                for i = 1:sfft
                    tperm[i] = _reverse_num(i - 1, power) + 1
                    rone[i] = round(
                        exp(one(T) * 2big(pi) * im * (i - 1) / sfft);
                        digits = prec + 16,
                        base = 2,
                    )
                    rone_conj[i] = conj(rone[i])
                end
            end
        end
        return new{T,numdims,dims}(size_fft, tab_permut, root_one, root_one_conj)
    end
end


function PrepareFftBig(size_fft, type::DataType; kwargs...)
    return PrepareFftBig(size_fft, one(type); kwargs...)
end

function PrepareFftBig(size_fft; kwargs...)
    return PrepareFftBig(size_fft, one(BigFloat); kwargs...)
end

function PrepareFftBig(s::Integer, x::T; kwargs...) where {T<:AbstractFloat}
    return PrepareFftBig((s,), x; kwargs...)
end

"""
$(SIGNATURES)
"""
function fftbig!(
    par::PrepareFftBig{T,NUMDIMS,DIMS},
    signal;
    flag_inv = false,
) where {T,NUMDIMS,DIMS}
    for (inddim, dim) in enumerate(DIMS)
        s = size(signal, dim)
        @assert par.size_fft[inddim] == s "size of dim $dim is $s while $(par.sizefft[inddim]) is waiting "
        s_div2 = div(s, 2)
        len = s
        n_len = len >> 1
        nb_r = 1
        rootO = flag_inv ? par.root_one[inddim] : par.root_one_conj[inddim]
        perm = par.tab_permut[inddim]
        while n_len != 0
            start = 1
            suite = start + n_len
            for i = 1:nb_r
                deb = 1
                for j = deb:nb_r:s_div2
                    if NUMDIMS == 1
                        signal[start], signal[suite] = (signal[start] + signal[suite]),
                        (signal[start] - signal[suite]) * rootO[j]
                    else
                        s_start = selectdim(signal, dim, start)
                        s_suite = selectdim(signal, dim, suite)
                        s_sum = s_start + s_suite
                        s_diff = (s_start - s_suite) * rootO[j]
                        s_start .= s_sum
                        s_suite .= s_diff
                    end
                    start += 1
                    suite += 1
                end
                start = suite
                suite = start + n_len
            end
            len = n_len
            n_len >>= 1
            nb_r <<= 1
        end
        signal .= if NUMDIMS == 1
            signal[perm]
        else
            selectdim(signal, dim, perm)
        end

        if flag_inv
            signal ./= s
        end
    end
    return signal
end

function fftbig(
    par::PrepareFftBig{T,NUMDIMS,DIMS},
    signal;
    flag_inv = false,
) where {T,NUMDIMS,DIMS}
    fl = flag_inv
    return fftbig!(
        par::PrepareFftBig,
        copy(convert(Array{Complex{T}}, signal));
        flag_inv = fl,
    )
end

fftgenall(_::Any, t::AbstractArray{Float64}) = fft(t, ntuple(x -> x, ndims(t)))

function fftgen!(
    _::PrepareFftBig{T,NUMDIMS,DIMS},
    t::AbstractArray{Complex{Float64}},
) where {T,NUMDIMS,DIMS}
    return fft!(t, DIMS)
end

function fftgen(
    _::PrepareFftBig{T,NUMDIMS,DIMS},
    t::AbstractArray{Float64},
) where {T,NUMDIMS,DIMS}
    return fft(t, DIMS)
end

fftgenall(p::PrepareFftBig, t::AbstractArray{Float64}) = fftgen(p, t)

fftgen(p::PrepareFftBig, t::AbstractArray{T}) where {T<:AbstractFloat} = fftbig(p, t)

function fftgen(p::PrepareFftBig, t::AbstractArray{Complex{T}}) where {T<:AbstractFloat}
    return fftbig(p, t)
end

function fftgen!(p::PrepareFftBig, t::AbstractArray{Complex{T}}) where {T<:AbstractFloat}
    return fftbig!(p, t)
end

fftgenall(p::PrepareFftBig, t::AbstractArray{T}) where {T} = fftgen(p, t)
fftgenall(p::PrepareFftBig, t::AbstractArray{Complex{T}}) where {T} = fftgen(p, t)
fftgenall!(p::PrepareFftBig, t::AbstractArray{Complex{T}}) where {T} = fftgen!(p, t)

ifftgenall(_::Any, t::AbstractArray{Complex{Float64}}) = ifft(t, ntuple(x -> x, ndims(t)))

function ifftgen(
    _::PrepareFftBig{T,NUMDIMS,DIMS},
    t::AbstractArray{Complex{Float64}},
) where {T,NUMDIMS,DIMS}
    return ifft(t, DIMS)
end

function ifftgen!(
    _::PrepareFftBig{T,NUMDIMS,DIMS},
    t::AbstractArray{Complex{Float64}},
) where {T,NUMDIMS,DIMS}
    return ifft!(t, DIMS)
end

ifftgenall(p::PrepareFftBig, t::AbstractArray{Complex{Float64}}) = ifftgen(p, t)

function ifftgen(p::PrepareFftBig, t::AbstractArray{Complex{T}}) where {T}
    return fftbig(p, t; flag_inv = true)
end

function ifftgen!(p::PrepareFftBig, t::AbstractArray{Complex{T}}) where {T}
    return fftbig!(p, t; flag_inv = true)
end

ifftgenall(p::PrepareFftBig, t::AbstractArray{Complex{T}}) where {T} = ifftgen(p, t)
ifftgenall!(p::PrepareFftBig, t::AbstractArray{Complex{T}}) where {T} = ifftgen!(p, t)
