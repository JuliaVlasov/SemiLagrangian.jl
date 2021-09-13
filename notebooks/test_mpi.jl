
using MPI
using DoubleFloats

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
            # for Float64 or Double64 ... for example
            for i = 1:mpid.nb
                vbcast = view(data, t_split[i])
                MPI.Bcast!(vbcast, i - 1, mpid.comm)
            end
        else
            # for BigFloat ... for example
            for i = 1:mpid.nb
                vbcast = view(data, t_split[i])
                bufr = MPI.bcast(vbcast, i - 1, mpid.comm)
                if i != mpid.ind
                    copy!(vbcast, bufr)
                end
            end
        end
     
    MPI.Barrier(mpid.comm)
end
function splititr(nb, lgtot)
    lg, r = divrem(lgtot, nb)
    return vcat(
        map(x -> ((x-1)*(lg+1)+1):x*(lg+1), 1:r),
        map(x -> ((x-1)*lg+r+1):(x*lg+r), (r+1):nb),
    )
end

function splitvec(nb, v)
    return map(x -> v[x], splititr(nb, length(v)))
end

function fct(T, sz::NTuple{2,Int})
    mpid = MPIData()
    nbsplit = mpid.nb

    t_itr = (splitvec(nbsplit, CartesianIndices(sz)),)

    @show typeof(t_itr)

    li = LinearIndices(sz)
    it = t_itr[1]
    tt_split = ntuple(y -> (li[it[y][1]]:li[it[y][end]]), nbsplit)

    @show typeof(tt_split)

    data = zeros(T,sz)
    cpt = 1
    loc_t = ref_t = time_ns()
    while true
        dec1 = rand()-0.5
        dec2 = rand()-0.5
        for i in tt_split[mpid.ind]
            a = div(i,sz[1])
            b = mod(i,sz[2])
            data[i] = (data[i] + dec1*a + dec2*b)/(a+b)
        end
        mpibroadcast(mpid, tt_split, data)
        if mpid.ind == 1
                if cpt % 80 == 0
                    n_t = time_ns()
                    loc = (n_t - loc_t)*1e-9
                    gen = (n_t - ref_t)*1e-9
                    moy = gen *80 / cpt
                    loc_t = n_t
                    println(". $(sum(data)) $cpt $moy $loc $gen")
                else
                    if cpt % 1 == 0
                        print(".")
                    end
                end
        end
        flush(stdout)
        cpt += 1
    end

end
function run_main(T,sz)
    fct(T,sz)   
end

run_main(Float64, (512,512) )
