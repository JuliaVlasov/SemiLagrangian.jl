






mutable struct RotationVar{T,N} <: AbstractExtDataAdv{T,N}
    decfl::Any
    decint::Any
    function RotationVar(adv::Advection{T,N}) where {T,N}
        return new{T,N}(missing,missing)
    end
end




"""
    initcoef!(pv::RotationVar{T, Nsp, Nv}, self::AdvectionData{T, Nsp, Nv, Nsum})

Implementation of the interface function that is called at the begining of each advection
    This is implementation for Vlasov-Poisson equation

"""

function initcoef!(
    pv::RotationVar{T,N},
    self::AdvectionData{T,N},
) where {T,N}
    st_cur, st_other = getst(self).perm
    mesh_cur = self.adv.t_mesh[st_cur]
    mesh_other = self.adv.t_mesh[st_other]
    sign = ( st_cur == 1 ) ? -1 : 1
    pv.decfl = sign * getcur_t(self)/step(mesh_cur) * mesh_other.points

    #  if isvelocitystate(self)
    #     #        println("sp trace init moins")
    #     pv.bufcur = (getcur_t(self) / step(mesh_v)) * mesh_sp.points
    # else
    #     mesh_sp = self.adv.t_mesh_sp[state_dim]
    #     #        println("sp trace init moins")
    #     pv.bufcur = (-getcur_t(self) / step(mesh_sp)) * mesh_v.points
    # end
end
function getrotationvar(adv::Advection)
    return RotationVar(adv)
end

"""
    getalpha(pv::RotationVar, self::AdvectionData, ind) 

Implementation of the interface function that is called before each interpolation in advection

"""
getalpha(pv::RotationVar, self::AdvectionData, ind) = (pv.decfl[ind],)
