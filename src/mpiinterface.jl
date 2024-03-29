using MPI
struct MPIData
    comm::Any
    nb::Any
    ind::Any
    function MPIData()
        if !MPI.Initialized()
            MPI.Init()
        end
        comm = MPI.COMM_WORLD
        nb = MPI.Comm_size(comm)
        ind = MPI.Comm_rank(comm) + 1
        return new(comm, nb, ind)
    end
end

function mpibroadcast(mpid, t_split, data::Array{T,N}) where {T,N}
    MPI.Barrier(mpid.comm)
    if T == Float64
        # for Float64 ... for example
        # Double64 here leads to memory leaks
        for i = 1:(mpid.nb)
            vbcast = view(data, t_split[i])
            MPI.Bcast!(vbcast, i - 1, mpid.comm)
        end
    else
        # for BigFloat ... for example
        # for Double64 also
        for i = 1:(mpid.nb)
            vbcast = view(data, t_split[i])
            bufr = MPI.bcast(vbcast, i - 1, mpid.comm)
            if i != mpid.ind
                copy!(vbcast, bufr)
            end
        end
    end
    return MPI.Barrier(mpid.comm)
end
