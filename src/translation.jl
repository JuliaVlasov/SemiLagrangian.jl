
mutable struct TranslationVar{T,N} <: AbstractExtDataAdv
    values::NTuple{N,T}
    valok::Any
end
function gettranslationvar(v::NTuple{N,T}) where {T,N}
    return TranslationVar{T,N}(v, ntuple(x -> zero(T), N))
end

function initcoef!(pv::TranslationVar{T,N}, self::AdvectionData{T,N}) where {T,N}
    st = getst(self)
    return pv.valok = ntuple(i -> pv.values[st.perm[i]] * getcur_t(self), st.ndims)
end

"""
    getalpha(pv::TranslationVar, self::AdvectionData, i, ind) 

Implementation of the interface function that is called before each interpolation in advection

"""
function getalpha(pv::TranslationVar{T,N}, self::AdvectionData{T,N}, ind) where {T,N}
    return pv.valok
end
