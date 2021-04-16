
function _clock()
    return ccall(:clock, Int32, ())
end
abstract type AbstractClockObs end
struct ClockObs{N} <: AbstractClockObs
    tbeg
    tsum
    function ClockObs(nb)
        tbeg = zeros(Int32, nb)
        tsum = zeros(Float64,nb)
        return new{nb}(tbeg,tsum)
    end
end

function clockbegin(cl::ClockObs, ind)
    cl.tbeg[ind] = _clock()
end
function clockend(cl::ClockObs, ind)
    cl.tsum[ind] += (_clock() - cl.tbeg[ind])/1000000
end
clocksum(cl::ClockObs, ind)=cl.tsum[ind]
struct NoClockObs <: AbstractClockObs
end
clockbegin(::AbstractClockObs, ind)=missing
clockend(::AbstractClockObs, ind)=missing
printall(::AbstractClockObs)=missing


function clockreset(cl::ClockObs)
    fill!(cl.tbeg,0)
    fill!(cl.tsum,0.0)
end

function printall(cl::ClockObs{N}) where {N}
    for i=1:N
        println("ClockObs i=$i t=$(cl.tsum[i])")
    end
end

    

