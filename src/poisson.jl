
include("fftbig.jl")
include("advection.jl")

"""

    compute_e!( e, mesh, ρ)

    ∇.e = - ρ

Inplace computation of electric field. Fields e and rho are
already allocated.

"""
function compute_e!(
    e::Vector{T},
    meshx::UniformMesh{T},
    rho::Vector{T},
) where {T}

    nx = meshx.length
    k = 2π / (meshx.stop - meshx.start)
    modes = zeros(Float64, nx)
    modes .= k .* vcat(0:nx÷2-1, -nx÷2:-1)
    modes[1] = 1.0
    e .= real(ifft(-1im .* fft(rho) ./ modes))

end
function compute_elfieldother(
    t_mesh_x::NTuple{N,UniformMesh{T}},
    rho::Array{T, N},
    pfft
) where {T <: AbstractFloat, N}

    fct_k(v)= im/sum(v.^2)

    v_k = vec_k_fft.(t_mesh_x)
    sz = length.(t_mesh_x)

    buf = fftgenall(pfft, rho) .* fct_k.(collect(Iterators.product(v_k...)))
    buf[1] = 0im
   
#    println("size(buf)=$(size(buf)) size(array_k)=$(size(array_k))")

    return ntuple( 
    x -> real(ifftgenall(pfft, reshape(v_k[x],tupleshape(x,N,sz[x])) .* buf )),
    N
)
end

function compute_elfield(
    t_mesh_x::NTuple{N,UniformMesh{T}},
    rho::Array{T, N},
    pfft
) where {T <: AbstractFloat, N}

    fct_k(ind,v)= v[ind] == 0 ? 0im : im*v[ind]/sum(v.^2)

    v_k = vec_k_fft.(t_mesh_x)

    buf = fftgenall(pfft, rho)

    println("buf[1]=$(buf[1])")
 
    array_k = collect(Iterators.product(v_k...))

#    println("size(buf)=$(size(buf)) size(array_k)=$(size(array_k))")

    return ntuple( 
        x -> real(ifftgenall(pfft, fct_k.(x, array_k) .* buf )),
        N
    )

end
# todo à optimiser
function compute_elfield!(
    e::NTuple{N,Array{Complex{T}, N}},
    t_mesh_x::NTuple{N,UniformMesh{T}},
    rho::Array{T, N},
    pfft
) where {T <: AbstractFloat, N}

    fct_k(ind,v)= v[ind] == 0 ? 0im : im*v[ind]/sum(v.^2)

    v_k = vec_k_fft.(t_mesh_x)

    buf = fftgenall(pfft, rho)

    println("buf[1]=$(buf[1])")
 
    array_k = collect(Iterators.product(v_k...))

#    println("size(buf)=$(size(buf)) size(array_k)=$(size(array_k))")

    for i=1:N
        e[i] .= fct_k.(i, array_k) .* buf
        ifftgenall!(pfft, e[i])
        e[i] .= real(e[i])
    end 

end





function compute_elfieldother!(
    t_e::NTuple{N,Array{Complex{T},N}},
    t_mesh_x::NTuple{N,UniformMesh{T}},
    rho::Array{T, N},
    pfft
) where {T <: AbstractFloat, N}
println("trace OK")
    res = compute_elfieldother(t_mesh_x,rho,pfft)
    for i=1:N
        t_e[i] .= res[i]
    end
    t_e
end
function compute_elfieldother2!(
    t_e::NTuple{N,Array{Complex{T},N}},
    t_mesh_x::NTuple{N,UniformMesh{T}},
    rho::Array{T, N},
    pfft
) where {T <: AbstractFloat, N}

    fct_k(v)= im/sum(v.^2)

    v_k = vec_k_fft.(t_mesh_x)
    sz = length.(t_mesh_x)

    buf = fftgenall(pfft, rho) .* fct_k.(collect(Iterators.product(v_k...)))
    buf[1] = 0im
   
#    println("size(buf)=$(size(buf)) size(array_k)=$(size(array_k))")
    for i=1:N
        t_e[i] .= reshape(v_k[i],tupleshape(i,N,sz[i])) .* buf
        t_e[i] .= real(ifftgenall(pfft,t_e[i]))
    end
    t_e
end

struct PoissonConst{T,Nsp,Nv}
    adv
    v_k
    fctv_k
    pfftbig
    function PoissonConst(
    adv::Advection{T, Nsp, Nv, Nsum}; 
    isfftbig=true
) where{T, Nsp, Nv, Nsum}
        Nsp <= Nv || thrown(ArgumentError("Nsp=$Nsp must less or equal to Nv=$Nv"))
        fct_k(v) = im / sum(v .^ 2)
        v_k = vec_k_fft.(adv.t_mesh_sp)
        sz = length.(adv.t_mesh_sp)
        fctv_k = fct_k.(collect(Iterators.product(v_k...)))
        pfftbig = if isfftbig
            PrepareFftBig(sz, T, numdims=N, ntuple(x->x,N))
        else
            missing
        end
        return new{T,Nsp,Nv}(adv, v_k, fctv_k, pfftbig)
    end
end

struct PoissonVar{T, Nsp, Nv}
    pc::PoissonConst{T, Nsp, Nv}
    rho::Array{T, Nsp}
    t_elfield::NTuple{Nsp,Array{Complex{T}, Nsp}}
    pfftbig
    bufcur
    function PoissonEq(
        pc::PoissonConst{T, Nsp, Nv}
) where{T, Nsp, Nv}
        sz = length.(pc.t_mesh_x)
        rho = Array{T, Nsp}(undef, sz)
        t_elfield = ntuple( x-> Array{T, N}(undef, sz), N)
        return new{T, Nsp, Nv}(rho, t_elfield, t_mesh_x, pfftbig, bufcur)
    end
end

function init!(self::AdvectionData{T, Nsp, Nv, Nsum}) where{T, Nsp, Nv, Nsum}
    pv::PoissonVar{T, Nsp, Nv} = getext(self)
    state_dim = getstate_dim(self)
    if state_dim == 1 && isvelocitystate(self)
        data = getdata(self)
        dv = prod(step, pv.pc.adv.t_mesh_v)
        dp.rho .= dv * reshape(sum(data, dims = ntuple(x->x+Nsp,Nv) ), size(dp.rho))
        dp.rho .-= sum(rho)/length(rho)

        buf = fftgenall(pfft, dp.rho) .* pv.pc.fctv_k
        buf[1] = 0im
       
    #    println("size(buf)=$(size(buf)) size(array_k)=$(size(array_k))")
    
        pv.t_elfield = ntuple( 
    x -> real(ifftgenall(pfft, reshape(pv.pc.v_k[x],tupleshape(x,Nsp,sz[x])) .* buf )),
    Nsp
)
    end
end

function getalpha(self::AdvectionData{T, Nsp, Nv, Nsum}, ind) where{T, Nsp, Nv, Nsum}
    pv::PoissonVar{T, Nsp, Nv} = getext(self)
    state_dim = getstate_dim(self)
    if isvelocitystate(self)
        return pv.t_elfield[state_dim][ind[1:Nsp]...]
    else
        mesh = pv.adv.t_mesh_v[state_dim]
        return (getcur_t(self) / mesh.step) * mesh.points[ind[state_dim+Nsp]]
    end
end
    

