using MPI

struct MPIData
    comm
    nb
    ind
    function MPIData()
        MPI.Init()
        comm = MPI.COMM_WORLD
        nb = MPI.Comm_size(comm)
        ind = MPI.Comm_rank(comm)+1
        return new(comm, nb, ind)
    end
end

function mpibroadcast(mpid, t_split, data::Array{T,N}) where {T,N}

    MPI.Barrier(mpid.comm)
    if isbitstype(T)
        # for Float64 or Double64 ... per example
        for i=1:mpid.nb
            vbcast = view(data, t_split[i])
            MPI.Bcast!(vbcast, i-1, mpid.comm)
        end
    else
        # for BigFloat ... per example
        for i=1:mpid.nb
            vbcast = view(data, t_split[i])
            bufr = MPI.bcast(vbcast, i-1, mpid.comm)
            if i != mpid.ind
                copy!(vbcast, bufr)
            end
        end
    end
    MPI.Barrier(mpid.comm)
end
