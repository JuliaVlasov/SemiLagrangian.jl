
include("../src/fctthreads.jl")

using Primes

nball=0

function fct(tab, ind, nb)
    m = ind-1
    for (i, n) in enumerate(tab)
        if i%nb == m       
            isprime(big"2"^n-1)
            global nball
            nball += 1
        end
    end
end



function test_fctthreads()
    nb = Threads.nthreads()-1

    data = zeros(Int64,nb)
    fill!(data,4423)

    ft = FctThreads(nb, fct, data)

    launchall(ft)

    for i = 1:10
        postall(ft)
        waitall(ft)
    end

    endall(ft)
end
function test_thrthreads()
    nb = Threads.nthreads()-1

    data = zeros(Int64,nb)
    fill!(data,4423)

    for i = 1:10
        Threads.@threads for ind = 1:nb
            fct(data, ind, nb)
        end
    end
end


println("nball=$nball")
@time test_fctthreads()
println("nball=$nball")
@time test_thrthreads()

println("nball=$nball")

