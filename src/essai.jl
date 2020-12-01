
include("advection.jl")
include("lagrange.jl")
include("interpolation.jl")
include("clockobs.jl")

cl_obs=ClockObs(5)

# modulo for "begin to one" array
modone(ind, n)=(n+ind-1)%n+1
function mkcoef(alpha, order)
    a = rand(order+1) .- 0.5
    mid = div(order+1,2)
    return map(i -> (i+alpha)*a[i+mid+1], -mid:order-mid)
end




function advbid!(bigbuf, order)
    T = Float64
    buf=Vector{T}(undef,size(bigbuf,1))
    coef = mkcoef(rand(), order)
    println("length(coef)=$(length(coef))")
    for i=1:size(bigbuf,2)
        f = view(bigbuf,:,i)
        interpolbid!(buf,f,coef)
        f .= buf
    end
end


function advbid2!(bigbuf, interp::InterpolationType{T,true}) where{T}
    order = get_order(interp)
    buf=Vector{T}(undef,size(bigbuf,1))
#    coef = mkcoef(rand(), order)
    alpha = rand()
    decint = convert(Int, floor(alpha))
    decfloat = alpha - decint
    precal = get_precal(interp,decfloat)
     println("length(precal)=$(length(precal)), typeof(precal)=$(typeof(precal))")
    for i=1:size(bigbuf,2)
        f = view(bigbuf,:,i)
        interpolbid2!(buf,f,decint,precal,interp)
        f .= buf
    end
end
function advbid!(bigbuf, interp::InterpolationType{T}) where{T}
    buf=Vector{T}(undef,size(bigbuf,1))
    alpha = rand()
    decint = convert(Int, floor(alpha))
    decfloat = alpha - decint
    precal = get_precal(interp, decfloat)
    println("length(precal)=$(length(precal))")

    for i=1:size(bigbuf,2)
        f = view(bigbuf,:,i)
        interpolate!(buf,f,decint,precal,interp)
        f .= buf
    end
end
function advbidth!(bigbuf, interp::InterpolationType{T}) where{T}
    buflgn=Array{T}(undef,size(bigbuf,1),4)
    alpha = rand()
    decint = convert(Int, floor(alpha))
    decfloat = alpha - decint
    precal = get_precal(interp, decfloat)
    println("length(precal)=$(length(precal))")
    sz = div(size(bigbuf,2),4)
    a = map(x->view(bigbuf,:,(x-1)*sz+1:x*sz), 1:4)
    Threads.@threads for ith=1:4
        buf = view(buflgn,:,ith)
        bb = a[ith]
        for i=1:sz
            f = view(bb,:,i)
            interpolate!(buf,f,decint,precal,interp)
            f .= buf
        end
    end
end
function advbidth2!(bigbuf, interp::InterpolationType{T}) where{T}
    buflgn=Array{T}(undef,size(bigbuf,1),4)
    alpha = rand()
    decint = convert(Int, floor(alpha))
    decfloat = alpha - decint
    precal = get_precal(interp, decfloat)
    println("length(precal)=$(length(precal))")
    Threads.@threads for i=1:size(bigbuf,2)
        buf = view(buflgn,:,Threads.threadid())
        f = view(bigbuf,:,i)
        interpolate!(buf,f,decint,precal,interp)
        f .= buf
    end
end
# function advbidth!(bigbuf, order)
#     T = Float64
#     buflgn=Array{T}(undef,size(bigbuf,1),4)
#     coef = mkcoef(rand(), order)
#     sz = div(size(bigbuf,2),4)

#     a = map(x->copy(view(bigbuf,:,(x-1)*sz+1:x*sz)), 1:4)

#     Threads.@threads for ith=1:4
#         println("threadid=$(Threads.threadid())")
#         buf = view(buflgn,:,ith)
#         bb = a[ith]
#         for i=1:sz
#             f = view(bb,:,i)
#             interpolbid!(buf,f,coef)
#             f .= buf
#         end
#     end
# end

T = BigFloat
@time bigbuf = rand(T, 32, 32^3)

@time for i= 1:1
    advbid!(bigbuf, 29)
end


interp = Lagrange(T, 29)
println("interp=$(get_type(interp))")


@time for i= 1:1
    advbid2!(bigbuf, interp)
end



@time for i= 1:1
    advbid!(bigbuf, interp)
end
@time for i= 1:1
    advbidth!(bigbuf, interp)
end
@time for i= 1:1
    advbidth2!(bigbuf, interp)
end

# printall(cl_obs)
# @time for i= 1:1
#     advbidth!(bigbuf, 29)
# end




