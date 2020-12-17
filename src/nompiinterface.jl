struct MPIData
    function MPIData()
        error("no MPI on this package")
    end
end

function mpibroadcast(mpid, t_split, data::Array{T,N}) where {T,N}
    error("no MPI on this package")
end

