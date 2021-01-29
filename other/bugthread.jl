
function fctbad(f, isthr)
    tabbuf = zeros(size(f,1),isthr ? Threads.nthreads() : 1)
    nb = size(f,2)
    if isthr
        @Threads.threads for i=1:nb
            buf = tabbuf[:,Threads.threadid()]
            buf .= 5f[:,i]
            f[:,i] .= buf
        end
    else
        buf = tabbuf[:,1]
        for i=1:nb
            buf .= 5f[:,i]
            f[:,i] .= buf
        end
    end
 
end
function fctgood(f, isthr)
    tabbuf = zeros(size(f,1),isthr ? Threads.nthreads() : 1)
    nb = size(f,2)
    if isthr
        @Threads.threads for i=1:nb
            buf = tabbuf[:,Threads.threadid()]
            buf .= 5f[:,i]
            f[:,i] .= buf
        end
    else
        buf2 = tabbuf[:,1]
        for i=1:nb
            buf2 .= 5f[:,i]
            f[:,i] .= buf2
        end
    end
 
end
data = zeros(2, 8)
for i=1:size(data,2), j=1:size(data,1)
    data[j,i] = i
end

println("avant data=$data")

fctbad(data, true)

println("apres data=$data")

data = zeros(2, 8)
for i=1:size(data,2), j=1:size(data,1)
    data[j,i] = i
end

println("avant data=$data")

fctgood(data, true)

println("apres data=$data")
