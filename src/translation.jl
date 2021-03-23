

struct TranslationVar{T,Nsum} <: AbstractExtDataAdv1d{T,Nsum}
    values::NTuple{Nsum,T}
end
function gettranslationvar(v::NTuple{Nsum,T}) where {T,Nsum}
    return TranslationVar{T,Nsum}(v)
end

"""
    getalpha(pv::TranslationVar, self::Advection1dData, ind) 

Implementation of the interface function that is called before each interpolation in advection

"""
function getalpha(
    pv::TranslationVar{T,Nsum},
    self::Advection1dData{T,Nsp,Nv,Nsum},
    ind,
) where {T,Nsp,Nv,Nsum}
    return pv.values[_getcurrentindice(self)] * getcur_t(self)
end
