
include("../src/fctthreads.jl")
include("../src/clockobs.jl")

using Primes

cl_obs=ClockObs(2)

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
    global cl_obs
    clockbegin(cl_obs,1)
    for i = 1:10
        postall(ft)
        waitall(ft)
    end
    clockend(cl_obs,1)

    endall(ft)
end
function test_thrthreads()
    nb = Threads.nthreads()

    data = zeros(Int64,nb)
    fill!(data,4423)

    global cl_obs
    clockbegin(cl_obs,2)
   for i = 1:10
        Threads.@threads for ind = 1:nb
            fct(data, ind, nb)
        end
    end
    global cl_obs
    clockend(cl_obs,2)

end


println("nball=$nball")
@time test_fctthreads()
println("nball=$nball")
@time test_thrthreads()

println("nball=$nball")
printall(cl_obs)

