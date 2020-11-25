



mutable struct FctThreads{N}
    sembegin
    semend
    fct
    data
    isend
    function FctThreads(nb, fct, data)
        sembegin = map( x-> Base.Semaphore(1), 1:nb)
        Base.acquire.(sembegin)
        semend = map( x-> Base.Semaphore(1), 1:nb)
        Base.acquire.(semend)
        isend = falses(nb)
        return new{nb}(sembegin, semend, fct, data, isend)
    end
end

function _launchonethread(ft::FctThreads{N}, ind) where{N}
    while true
        Base.acquire(ft.sembegin[ind])
        if ft.isend[ind]
            break
        end
        ft.fct(ft, ind, N)
        Base.release(ft.semend[ind])
    end
end

function launchall(ft::FctThreads{N}) where{N}
    for i=1:N
        Threads.@spawn _launchonethread(ft, i)
    end
end
function postall(ft::FctThreads{N}) where{N}
    for i=1:N
        Base.release(ft.sembegin[i])
    end
end
function waitall(ft::FctThreads{N}, isend=false) where{N}
    for i=1:N
        Base.acquire(ft.semend[i])
        if isend
            ft.isend[i] = true
            Base.release(ft.sembegin[i])
        end
    end
end
function endall(ft::FctThreads{N}) where{N}
    for i=1:N
        ft.isend[i] = true
        Base.release(ft.sembegin[i])
    end
end



