using MPI
MPI.Init()

function getview(data, id, n)
    @assert mod(length(data),n) == 0 "length(data)=$(length(data)) must be null modulus n=$n"
    s = div(length(data),n)
    deb=id*s+1
    fin=deb+s-1
    return view(data,deb:fin)
end

function myfctalltoall(data,comm)
    n = MPI.Comm_size(comm)
    id = MPI.Comm_rank(comm)
    vsend = getview(data, id, n)

    for i=0:n-1
        if i != id
            MPI.Send(vsend, i, 9999, comm)
        end
    end
    for i=0:n-1
        if i != id
            vrecv = getview(data, i, n)
            MPI.Recv!(vrecv, i, 9999, comm)
        end
    end
    MPI.Barrier(comm)
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
    T = Float64
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
end



comm = MPI.COMM_WORLD

T = BigFloat
bigbuf = rand(T, 64, 32^3 * 3)

@time for i= 1:1
    advbid!(bigbuf, 29, comm)
end

# n = MPI.Comm_size(comm)
# id = MPI.Comm_rank(comm)
# println("Hello world, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")
# s = 4
# A = zeros( n*s)
# deb=id*s+1
# fin=deb+s-1
# for i=deb:fin
#         A[i] = 100*id+i
# end

# println("avant alltoall I am $id A=$A")

# myfctalltoall(A,comm)

# println("apres alltoall I am $id A=$A")

