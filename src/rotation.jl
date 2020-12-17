
include("advection.jl")

function _get_perm(adv::Advection{T, Nsp, Nv, Nsum, timeopt}, curstate) where {T, Nsp, Nv, Nsum, timeopt}
    return if isvelocity(adv, curstate)
        [2, 1]
    else
        [1, 2]
    end

end
function _get_split(adv, curstate)
    if adv.nbsplit != 1
        perm = _get_perm(adv,curstate)
        return splititr(adv.nbsplit, sizeall(adv)[perm][end])
    else
        return missing
    end
end
function _get_t_itrfirst(adv::Advection{T,Nsp, Nv, Nsum, timeopt}, t_itr, curid) where{T,Nsp,Nv,Nsum,timeopt}
    
    szitr = sizeitr(adv)
    if adv.nbsplit == 1
        if isvelocity(adv, curid)
#            println("_get_t_itrfirst trace1")
            return Iterators.product(szitr[1:Nsp]...)
        else
#            println("_get_t_itrfirst trace2 return=$(szitr[curid+Nsp])")
            return szitr[curid+Nsp]
        end
    else
        if isvelocity(adv, curid)
#            println("_get_t_itrfirst trace3")
            return ntuple( x -> Iterators.product((szitr[1:Nsp-1]..., (t_itr[x],)...)...), adv.nbsplit)
        else
#            println("_get_t_itrfirst trace4")
            return t_itr
        end
    end
end



struct RotationConst{T, Nsp, Nv}
    adv
    t_perms
    tt_split
    t_itrfirst
    t_itrsecond
    function RotationConst(
    adv::Advection{T, Nsp, Nv, Nsum, timeopt}
) where{T, Nsp, Nv, Nsum, timeopt}
        Nsp == Nv || thrown(ArgumentError("Nsp=$Nsp must be equal to Nv=$Nv"))
        t_perms = ntuple(x -> _get_perm(adv, x), Nsum)
        tt_split = ntuple(x-> _get_split(adv, x), Nsum)
        t_itrfirst = ntuple(x -> _get_t_itrfirst(adv, tt_split[x], x), Nsum)
        t_itrsecond = ntuple(x -> Iterators.product(sizeitr(adv)[t_perms[x]][2:Nsum-1]...), Nsum)
        return new{T,Nsp,Nv}(adv, t_perms, tt_split, t_itrfirst, t_itrsecond)
    end
end
getperm(pc::RotationConst,advd::AdvectionData)=pc.t_perms[_getcurrentindice(advd)]
gett_split(pc::RotationConst, advd::AdvectionData)=pc.tt_split[_getcurrentindice(advd)]


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
gett_split(pvar::RotationVar,advd::AdvectionData)=gett_split(pvar.pc,advd)



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







#Obsolete now
"""
    getalpha(self::AdvectionData{T, Nsp, Nv, Nsum}, ind) 

Implementation of the interface function that is called before each interpolation in advection

"""
function getalpha(self::AdvectionData{T, Nsp, Nv, Nsum}, ind, indthread) where{T, Nsp, Nv, Nsum}
    pv::PoissonVar{T, Nsp, Nv} = getext(self)
    state_dim = getstate_dim(self)
    alpha = if isvelocitystate(self)
        if indthread <= 0
            pv.bufcur[ind[1:Nsp]...]
        else
            pv.bufcur[ind[1:(Nsp-1)]...,((indthread-1)*self.adv.szsplit + ind[Nsp],)...]
        end
    else
        pv.bufcur[ind[self.state_dim+Nsp]]
    end
#    println("ind=$ind alpha=$alpha")
    return alpha
end
getalpha(self::AdvectionData, ind)=getalpha(self,ind,0)

function getprecal(pv::RotationVar{T, Nsp, Nv}, self::AdvectionData{T, Nsp, Nv, Nsum}, ind) where {T, Nsp, Nv, Nsum}
#    alpha = isvelocitystate(self) ? pv.bufcur[ind...] : pv.bufcur[ind]
    alpha = pv.bufcur[ind...]
    decint = convert(Int, floor(alpha))
    decfloat = alpha - decint
    return decint, get_precal(getinterp(self),decfloat)
end


function getitrfirst(pc::RotationConst, advd::AdvectionData{T,Nsp, Nv, Nsum, timeopt}) where{T,Nsp,Nv,Nsum,timeopt}
    itrfirst = pc.t_itrfirst[_getcurrentindice(advd)]
    if pc.adv.nbsplit != 1
        ind = timeopt == MPIOpt ? advd.adv.mpid.ind : Threads.threadid()
        return itrfirst[ind]
    else
#        println("trace good itrfirst=$itrfirst")
        return itrfirst
    end
 
end
getitrfirst(pvar::RotationVar, advd::AdvectionData)=getitrfirst(pvar.pc, advd)
addcolindend(ind::Tuple,tup)=(:,tup...,ind...)
addcolindend(ind::Int,tup)=(:,tup...,ind)

function getitrsecond(pc::RotationConst, advd::AdvectionData{T,Nsp, Nv, Nsum}, indfirst) where{T,Nsp,Nv,Nsum}
    perm = getperm(pc,advd)
    szitr = sizeitr(advd.adv)[perm]
    tupmid = isvelocitystate(advd) ? szitr[2:Nv] : szitr[2:Nsum-1]
    return addcolindend.((indfirst,), Iterators.product(tupmid...))
end
getitrsecond(pvar::RotationVar, advd::AdvectionData, indfirst)=getitrsecond(pvar.pc, advd, indfirst)
    
