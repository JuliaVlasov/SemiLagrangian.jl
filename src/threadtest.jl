global _nbtot=0
global tch = Vector{Channel}(undef,32)
global tab =zeros(Int64,1000,1000)


function fctlgn!(v,j)
    for i=1:length(v)
        v[i] += i +1+v[j%size(v,1)+1]
        v[i] = v[i] ^2
        global _nbtot
        _nbtot += 1
        v[i] %= 1021
    end
end


function fct!(ich)

    global tab
    s = div(size(tab,2),nb)
    local deb = (ich-1)*s+1
    local fin = deb+s-1

    println("fct2 itr=$itr")
    for i in deb:fin
        for j = 1:10000
            fctlgn!(tab[:,i],j)
        end
    end
    global tch
    put!(tch[ich],"fin")
end




function fctmain2(nb)

    tab = zeros(Int,100,1000)
    s = div(size(tab,2),nb)
    tch = Vector{Channel}(undef,nb)
    c_fct = @cfunction(fct!, Int32, (Ref{Int32}, ));
    for i=1:nb
        if nb != 1
#            Base.Threads.@spawn fct2!(tab, deb:fin, tch[i])
#            Base.Threads.@spawn fct!(par)
            tch[i] = Channel(1)
            pthr = Vector{Int32}(undef,4)
            println("trace98 i=$i")
            ccall(:pthread_create, Int32, (Ptr{Int32}, Ptr{Int32}, Ptr{Cvoid}, Ref{Int32}), pthr, C_NULL, c_fct, i)
        else
#            fct2!(tab, deb:fin,tch[1])
            fct!(par)
        end
    end
    for i = 1:nb
        wait(tch[i])
    end

end

@time fctmain2(2)

println("$_nbtot")


