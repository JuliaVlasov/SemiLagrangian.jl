struct MPIData
    function MPIData()
        return error("no MPI on this package")
    end
end

function mpibroadcast(mpid, t_split, data::Array{T,N}) where {T,N}
    return error("no MPI on this package")
end
