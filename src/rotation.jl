


function _get_perm(adv::Advection{T, Nsp, Nv, Nsum, timeopt}, curstate) where {T, Nsp, Nv, Nsum, timeopt}
    return if isvelocity(adv, curstate)
        [2, 1]
    else
        [1, 2]
    end

end


struct RotationConst{T, Nsp, Nv}
    adv
    t_perms
    function RotationConst(
    adv::Advection{T, Nsp, Nv, Nsum, timeopt}
) where{T, Nsp, Nv, Nsum, timeopt}
        Nsp == Nv || thrown(ArgumentError("Nsp=$Nsp must be equal to Nv=$Nv"))
        t_perms = ntuple(x -> _get_perm(adv, x), Nsum)
        return new{T,Nsp,Nv}(adv, t_perms)
    end
end
getperm(pc::RotationConst,advd::AdvectionData)=pc.t_perms[_getcurrentindice(advd)]
getperm(pc::RotationConst,curst::Int)=pc.t_perms[curst]


mutable struct RotationVar{T, Nsp, Nv}
    pc::RotationConst{T, Nsp, Nv}
    bufcur
    function RotationVar(pc::RotationConst{T, Nsp, Nv}) where{T, Nsp, Nv}
        sz = length.(pc.adv.t_mesh_sp)
        rho = Array{T, Nsp}(undef, sz)
        return new{T, Nsp, Nv}(pc, missing)
    end
end
getperm(pvar::RotationVar,advd::AdvectionData)=getperm(pvar.pc,advd)
getperm(pvar::RotationVar,curst::Int)=getperm(pvar.pc, curst)




"""
    initcoef!(pv::RotationVar{T, Nsp, Nv}, self::AdvectionData{T, Nsp, Nv, Nsum})

Implementation of the interface function that is called at the begining of each advection
    This is implementation for Vlasov-Poisson equation

"""

function initcoef!(pv::RotationVar{T, Nsp, Nv}, self::AdvectionData{T, Nsp, Nv, Nsum}) where{T, Nsp, Nv, Nsum}
    state_dim = getstate_dim(self)
    mesh_sp = self.adv.t_mesh_sp[state_dim]
    mesh_v = self.adv.t_mesh_v[state_dim]
    if isvelocitystate(self)
#        println("sp trace init moins")
        pv.bufcur = (getcur_t(self)/step(mesh_v))*mesh_sp.points
    else
        mesh_sp = self.adv.t_mesh_sp[state_dim]
#        println("sp trace init moins")
       pv.bufcur = (-getcur_t(self)/step(mesh_sp))*mesh_v.points
    end
end
function getrotationvar(adv::Advection)
    pc = RotationConst(adv)
    return RotationVar(pc)
end

"""
    getalpha(self::AdvectionData{T, Nsp, Nv, Nsum}, ind) 

Implementation of the interface function that is called before each interpolation in advection

"""
getalpha(pv::RotationVar, self::AdvectionData, ind)=pv.bufcur[ind]


