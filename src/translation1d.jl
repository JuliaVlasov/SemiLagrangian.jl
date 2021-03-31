struct TranslationVar1d{T,Nsum} <: AbstractExtDataAdv1d{T,Nsum}
    values::NTuple{Nsum,T}
end
function gettranslationvar1d(v::NTuple{Nsum,T}) where {T,Nsum}
    return TranslationVar1d{T,Nsum}(v)
end
"""
    getalpha(pv::TranslationVar1d, self::Advection1dData, ind) 

Implementation of the interface function that is called before each interpolation in advection

"""
function getalpha(
    pv::TranslationVar1d{T,Nsum},
    self::Advection1dData{T,Nsp,Nv,Nsum},
    ind,
) where {T,Nsp,Nv,Nsum}
    return pv.values[_getcurrentindice(self)] * getcur_t(self)
end
