"""
$(TYPEDEF)
"""
mutable struct RotationVar{T,N} <: AbstractExtDataAdv
    decfl::Any
    decint::Any
    function RotationVar(adv::Advection{T,N}) where {T,N}
        return new{T,N}(missing, missing)
    end
end

"""
$(SIGNATURES)

    initcoef!(pv::RotationVar{T, Nsp, Nv}, self::AdvectionData{T, Nsp, Nv, Nsum})

Implementation of the interface function that is called at the begining of each advection
    This is implementation for Vlasov-Poisson equation

"""
function initcoef!(
    pv::RotationVar{T,N},
    self::AdvectionData{T,N,timeopt,NoTimeAlg},
) where {T,N,timeopt}
    st_cur, st_other = getst(self).perm
    mesh_cur = self.adv.t_mesh[st_cur]
    mesh_other = self.adv.t_mesh[st_other]
    sign = (st_cur == 1) ? -1 : 1
    return pv.decfl = sign * getcur_t(self) / step(mesh_cur) * mesh_other.points

end

"""
$(SIGNATURES)
"""
function initcoef!(
    pv::RotationVar{T,2},
    self::Union{
        AdvectionData{T,2,timeopt,ABTimeAlg_ip},
        AdvectionData{T,2,timeopt,ABTimeAlg_new},
    },
) where {T,timeopt}
    adv = self.adv
    sz = sizeall(adv)
    buf1 = -getcur_t(self) / step(adv.t_mesh[1]) * adv.t_mesh[2].points
    buf2 = getcur_t(self) / step(adv.t_mesh[2]) * adv.t_mesh[1].points
    if ismissing(self.bufcur)
        self.bufcur = zeros(OpTuple{2,T}, sz)
    end
    for ind in CartesianIndices(sizeall(adv))
        self.bufcur[ind] = OpTuple((buf1[ind.I[2]], buf2[ind.I[1]]))
    end

end

"""
$(SIGNATURES)
"""
function getrotationvar(adv::Advection)
    return RotationVar(adv)
end

"""
$(SIGNATURES)

    getalpha(pv::RotationVar, self::AdvectionData, ind) 

Implementation of the interface function that is called before each interpolation in advection

"""
getalpha(pv::RotationVar, self::AdvectionData, ind) = (pv.decfl[ind],)
