
function _clock()
    return ccall(:clock, Int32, ())
end
function _time()
    vtime = Vector{Int32}(undef,4)
    ccall(:gettimeofday, Int32, (Ptr{Int32}, Ptr{Int32}), vtime, C_NULL)
    return Int64(vtime[1])*1000000 + Int64(vtime[3])
end
struct ClockObs{N}
    cl_beg
    cl_sum
    t_beg
    t_sum
    function ClockObs(nb)
        cl_beg = zeros(Int32, nb)
        cl_sum = zeros(Float64,nb)
        t_beg = zeros(Int64, nb)
        t_sum = zeros(Float64,nb)
        return new{nb}(cl_beg, cl_sum, t_beg, t_sum)
    end
end
function clockbegin(cl::ClockObs, ind)
    cl.cl_beg[ind] = _clock()
    cl.t_beg[ind] = _time()
end
function clockend(cl::ClockObs, ind)
    cl.cl_sum[ind] += (_clock() - cl.cl_beg[ind])/1000000
    cl.t_sum[ind] += (_time() - cl.t_beg[ind])/1000000
end
clocksum(cl::ClockObs, ind)=cl.cl_sum[ind]
timesum(cl::ClockObs, ind)=cl.t_sum[ind]
function clockreset(cl::ClockObs)
    fill!(cl.cl_beg,0)
    fill!(cl.cl_sum,0.0)
    fill!(cl.t_beg,0)
    fill!(cl.t_sum,0.0)
end

function printall(cl::ClockObs{N}) where {N}
    for i=1:N
        println("ClockObs i=$i cl=$(cl.cl_sum[i]) t=$(cl.t_sum[i])")
    end
end

    

