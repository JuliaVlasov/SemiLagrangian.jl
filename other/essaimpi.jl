using MPI

MPI.Init()
using Polynomials
using DoubleFloats

function getview(data, id, n)
    @assert mod(length(data),n) == 0 "length(data)=$(length(data)) must be null modulus n=$n"
    s = div(length(data),n)
    deb=id*s+1
    fin=deb+s-1
    return view(data,deb:fin)
end
function getview(data, id, n, dim)
    @assert mod(size(data,dim),n) == 0 "length(data)=$(length(data)) must be null modulus n=$n"
    s = div(size(data,dim),n)
    deb=id*s+1
    fin=deb+s-1
    vd = view(data,1:length(data))
    return selectdim(vd, dim, deb:fin)
end

function myfctalltoallold(data,comm)
    n = MPI.Comm_size(comm)
    id = MPI.Comm_rank(comm)
    vsend = getview(data, id, n)

    vr = Vector{MPI.Request}(undef,(n-1)*2)

    ireq  =1

    for i=0:n-1
        if i != id
            vrecv = getview(data, i, n)
            vr[ireq] = MPI.Irecv!(vrecv, i, 9999, comm)
            ireq += 1
        end
    end

    println("trace3 from $id")

    for i=0:n-1
        if i != id
            vr[ireq] = MPI.Isend(vsend, i, 9999, comm)
            ireq += 1
        end
    end

    println("trace4 from $id")

  
    stats = MPI.Waitall!(vr)
    println("trace5 from $id")

    MPI.Barrier(comm)
    println("trace6 from $id")

end
function myfctalltoall(data::Array{Double64},comm)
    n = MPI.Comm_size(comm)
    MPI.Barrier(comm)
    for i=0:n-1
        vbcast = getview(data, i, n)
        MPI.Bcast!(vbcast, i, comm)
    end
end
function myfctalltoall(data::Array{Float64},comm)
    n = MPI.Comm_size(comm)
    MPI.Barrier(comm)
    for i=0:n-1
        vbcast = getview(data, i, n)
        MPI.Bcast!(vbcast, i, comm)
    end
end
function myfctalltoall(data::Array{BigFloat},comm)
    n = MPI.Comm_size(comm)
    id = MPI.Comm_rank
    MPI.Barrier(comm)
    s = div(length(data),n)
    for i=0:n-1
        vbcast = getview(data, i, n)
        bufr = MPI.bcast(vbcast, i, comm)
        if i != id
            copy!(vbcast, bufr)
        end
    end
end

modone(ind, n)=(n+ind-1)%n+1
function mkcoef(alpha, order)
    a = rand(order+1) .- 0.5
    mid = div(order+1,2)
    return map(i -> (i+alpha)*a[i+mid+1], -mid:order-mid)
end
function interpolbid!(lgnout, lgnin, coef)
    lg = length(lgnout)
    order = length(coef)-1
    mid = div(order+1,2)
    for i=1:lg
        indbeg = i-mid
 #       indend = i+order-mid
        indend = indbeg+order
        lgnout[i] = sum(lgnin[modone.(indbeg:indend,lg)] .* coef)
    end
end


function advbid!(bigbuf, order, comm)
    T = typeof(bigbuf[1])
    n = MPI.Comm_size(comm)
    id = MPI.Comm_rank(comm)
    szlgn = size(bigbuf,1)
    buf=Vector{T}(undef,szlgn)
    coef = mkcoef(rand(), order)
    println("length(coef)=$(length(coef))")
    myview = getview(bigbuf, id, n)
    nblgn = div(length(myview),szlgn)
    myview = reshape(myview, szlgn, nblgn)
    println("nblgn=$nblgn")
    for i=1:nblgn
        f = view(myview,:,i)
        interpolbid!(buf,f,coef)
        f .= buf
    end
    println("trace1 from $id")
    myfctalltoall(bigbuf,comm)
    println("trace2 from $id")
end



@time comm = MPI.COMM_WORLD

T = Double64
@time bigbuf = rand(T, 32, 32, 32, 16*3)

@time for i= 1:1
    advbid!(bigbuf, 29, comm)
end
# T = BigFloat
# comm = MPI.COMM_WORLD
# n = MPI.Comm_size(comm)
# id = MPI.Comm_rank(comm)
# println("Hello world, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")
# s = 4
# A = zeros(T, n*s)
# deb=id*s+1
# fin=deb+s-1
# for i=deb:fin
#         A[i] = 100*id+i
# end

# println("avant alltoall I am $id A=$A")

# myfctalltoall(A,comm)

# println("apres alltoall I am $id A=$A")

