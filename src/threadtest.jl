global _nbtot=0
global channelnum= Channel{Int}(1)
put!(channelnum,1)
global tch = Vector{Channel}(undef,32)
global tabthr = Vector{Int32}(undef,4*32)
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


function fct!(voidpar)
println("trace fct!1")
    global tab
    global channelnum
    ich = take!(channelnum)
    put!(channelnum,ich+1)

    s = div(size(tab,2),nb)
    local deb = (ich-1)*s+1
    local fin = deb+s-1

    println("fct2 itr=$ich")
    for i in deb:fin
        for j = 1:10000
            fctlgn!(tab[:,i],j)
        end
    end
    global tch
    put!(tch[ich],"fin")
end




function fctmain2(nb)
    global tab
    global tabthr
    s = div(size(tab,2),nb)
    tch = Vector{Channel}(undef,nb)
    c_fct = @cfunction(fct!, Any, (Ptr{Cvoid}, ));
    println("nb=$nb")
    for i=1:nb
        if nb != 1
#            Base.Threads.@spawn fct2!(tab, deb:fin, tch[i])
#            Base.Threads.@spawn fct!(par)
            tch[i] = Channel(1)
            indthr = (i-1)*4+1
            pthr = pointer(tabthr[indthr:indthr+3])
            println("trace98 i=$i")
            ccall(:pthread_create, Int32, (Ptr{Int32}, Ptr{Int32}, Ptr{Cvoid}, Ptr{Cvoid}), pthr, C_NULL, c_fct, pointer(tab))
        else
#            fct2!(tab, deb:fin,tch[1])
            fct!(par)
        end
    end
    for i = 1:nb
        wait(tch[i])
    end

end

@time fctmain2(30)

println("$_nbtot")


