
struct GeoConst{T,N}
    adv::Advection{T,N}
    #    rho::T
    #    g::T
    #    buoyancy_freq_n::T
    odg_b::T
    #    hv_order::Int
    #    hv_val::T
    coefrsqk::NTuple{N,Array{Complex{T},N}}
    pfftbig::Union{Missing,PrepareFftBig{T}}
    function GeoConst(
        adv::Advection{T,N};
        #        rho::T = 1e3,
        #        g::T = 9.81,
        buoyancy_freq_n::T = 3 * 2 * T(2pi) / (24 * 3600) * sin(T(pi) / 4),
        odg_b::T = T(1e-3),
        #        hv_order::Int = 8,
        #        hv_val::T = T(0),
    ) where {T,N}
        @assert N == 2 "the number of dimension must be 2"
        sz = sizeall(adv)
        coefrsqk = ntuple(x -> zeros(Complex{T}, sz), N)
        pfftbig = if T != Float64
            PrepareFftBig(sizeall(adv), T; numdims = N, dims = ntuple(x -> x, N))
        else
            missing
        end
        veck = vec_k_fft.(adv.t_mesh)
        for ind in CartesianIndices(sz)
            tp = ntuple(i -> veck[i][ind.I[i]], N)
            k = sqrt(sum(tp .^ 2))
            v = k == 0 ? 0 : 1 / k
            for d = 1:N
                coefrsqk[d][ind] = v .* tp[d%2+1]
            end
        end
        coefrsqk[1] .*= im # /buoyancy_freq_n
        coefrsqk[2] .*= -im # /buoyancy_freq_n

        return new{T,N}(
            adv,
            #            rho,
            #            g,
            #            buoyancy_freq_n,
            odg_b,
            #            hv_order,
            #            hv_val,
            coefrsqk,
            pfftbig,
        )
    end
end

struct GeoVar{T,N} <: AbstractExtDataAdv
    gc::GeoConst{T,N}
end

function getgeovar(adv::Advection{T,N}; kwargs...) where {T,N}
    return GeoVar{T,N}(GeoConst(adv; kwargs...))
end

function initdata!(geoc::GeoVar{T,N}, advd::AdvectionData{T,N}) where {T,N}
    mesh_x = advd.adv.t_mesh[1]
    mesh_y = advd.adv.t_mesh[2]
    lx = stop(mesh_x) - start(mesh_x)
    ly = stop(mesh_y) - start(mesh_y)
    x, y = points.(advd.adv.t_mesh)
    ee = 4
    σ = lx / 15 # Length scale close to the Rossby radius

    ## Spatial buoyancy field
    function anticyclone(cx, cy)
        return (
            exp.(-ee .* (x .- cx) .^ 2 ./ 2σ^2) .* transpose(exp.(-(y .- cy) .^ 2 ./ 2σ^2))
        )
    end

    # Warm anticyclones
    advd.data = anticyclone(lx / 4, ly / 4)
    advd.data .+= anticyclone(3lx / 4, ly / 4)

    # Cold cyclones
    advd.data .-= anticyclone(lx / 4, 3ly / 4)
    advd.data .-= anticyclone(3lx / 4, 3ly / 4)

    # Specify the amplitude of the buoyancy
    return advd.data .*= geoc.gc.odg_b
end

function initcoef!(geoc::GeoVar{T,N}, advd::AdvectionData{T,N}) where {T,N}
    pfft = geoc.gc.pfftbig
    dt = getcur_t(advd)
    buf = fftgenall(pfft, advd.data)

    if ismissing(advd.bufcur)
        advd.bufcur = zeros(OpTuple{N,T}, sizeall(advd.adv))
    end

    return advd.bufcur .=
        dt *
        OpTuple.(zip(ntuple(x -> real(ifftgenall(pfft, geoc.gc.coefrsqk[x] .* buf)), N)...))
end
@inline function getalpha(
    ::GeoVar{T,N},
    self::AdvectionData{T},
    indext::CartesianIndex,
    indbuf::CartesianIndex,
) where {T,N,Nsp,Nv}
    st = getst(self)
    ind = CartesianIndex((indbuf.I..., indext.I...)[st.invp])
    return ntuple(x -> self.bufcur[ind][st.invp[x]], st.ndims)
end
