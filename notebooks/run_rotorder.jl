include("../src/mesh.jl")
include("../src/advection.jl")
include("../src/lagrange.jl")
# include("../src/lagrange2d.jl")
# include("../src/advection2d.jl")
include("../src/bspline.jl")
include("../src/bsplinelu.jl")
include("../src/bsplinefft.jl")
include("../src/interpolation.jl")

using LinearAlgebra
using Plots


function exact(tf::T, mesh1::UniformMesh{T}, mesh2::UniformMesh{T}) where {T}

    f = zeros(T,(mesh1.length,mesh2.length))
    for (i, x) in enumerate(mesh1.points), (j, y) in enumerate(mesh2.points)
        xn = cos(tf) * x + sin(tf) * y
        yn = - sin(tf) * x + cos(tf) * y
        f[i,j] = exp(-6*((xn)^2+(yn+big"1.8")^2))
    end

    return f

end

function fctadv( interp, mesh1, mesh2, v1, v2, refdeb, refend, dt)
    adv_x1 = Advection( mesh1, interp )
    adv_x2 = Advection( mesh2, interp )

    f = copy(refdeb)
    ft = copy(refdeb)
    advection!(adv_x1, f,  v1, tan(dt/2))
    transpose!(ft, f)
    advection!(adv_x2, ft, v2, sin(dt))
    transpose!(f, ft)
    advection!(adv_x1, f, v1, tan(dt/2))

    return norm(f-refend,Inf)
end




function fctmain( sz, dt::T, ordmax) where{T}

    mesh1 = UniformMesh(T(-5), T(5), sz; endpoint=false)
    mesh2 = UniformMesh(T(-5), T(5), sz; endpoint=false)

    refdeb = exact( zero(T), mesh1, mesh2)
    refend = exact( dt, mesh1, mesh2)

    v1 = - collect(mesh2.points)
    v2 = + collect(mesh1.points)

    indmax=div(ordmax-3,2)+1

    y = ones(indmax,3)
    x = zeros(Int64,indmax)

    labels = Array{String, 2}(undef,1,3)
    labels[1,1] = "Lagrange"
    labels[1,2] = "B-Spline LU"
    labels[1,3] = "B-Spline FFT"

    ind=1

    for order=3:2:ordmax
        println("order=$order")
        lag = Lagrange(BigFloat, order)
        splu = B_SplineLU(order, mesh1.length, zero(T))
        spfft = B_SplineFFT(order, mesh1.length, zero(T))

        y[ind,1] = fctadv(lag, mesh1, mesh2, v1, v2, refdeb, refend, dt)
        y[ind,2] = fctadv(splu, mesh1, mesh2, v1, v2, refdeb, refend, dt)
        y[ind,3] = fctadv(spfft, mesh1, mesh2, v1, v2, refdeb, refend, dt)

        println("ind=$ind order=$order y[ind,:]=$(y[ind,:])")

        x[ind] = order


        ind += 1
    end

    p = Plots.plot(
        x,
        log10.(y),
        xlabel="order",
        ylabel="error",
        legend=:bottomleft,
        label=labels,
        marker=2
    )
    prec = precision(BigFloat)
    
    Plots.savefig(p, "out/result2_$(sz)_$(prec)_$(ordmax).pdf")
end

function fctinter( interp::InterpolationType, sz)

    fct(x)=cos(2big(pi)*x)
    mesh = big.(1:sz)/sz
    deb = fct.(mesh)
    dec = big"0.351726155665655665187291927162514231451"
    ref = fct.( mesh .+ dec/sz)
    res = zeros(BigFloat, sz)
    interpolate!(res, deb, dec, interp)
    return Float64(norm(ref-res,Inf))
end


function fctmain2(sz, ordmax)

    indmax=ordmax-2

    y = ones(indmax,2)
    x = zeros(Int64,indmax)

    labels = Array{String, 2}(undef,1,2)
    labels[1,1] = "Lagrange"
    labels[1,2] = "B-Spline LU"
#    labels[1,3] = "B-Spline FFT"

    ind=1

    for order=3:ordmax
        println("order=$order")
        lag = Lagrange(BigFloat, order)
        splu = B_SplineLU(order, sz, big"0.")
#        spfft = B_SplineFFT(order, sz, big"0.")

        y[ind,1] = fctinter(lag, sz)
        y[ind,2] = fctinter(splu, sz)
#        y[ind,3] = fctinter(spfft, sz)

        println("ind=$ind order=$order y[ind,:]=$(y[ind,:])")

        x[ind] = order


        ind += 1
    end

    p = Plots.plot(
        x,
        log10.(y),
        xlabel="order",
        ylabel="error",
        legend=:bottomleft,
        label=labels,
        marker=2
    )
    prec = precision(BigFloat)
    
    Plots.savefig(p, "out/resinterp_$(sz)_$(prec)_$(ordmax).pdf")
end

function fctmain_gnuplot(sz, ordmax)

    prec = precision(BigFloat)
    println("# begin gnuplot size=$sz julia precision=$prec")
    println("# order\tLagrange\tB_SplineLU")
    for order=3:ordmax
        lag = Lagrange(BigFloat, order)
        splu = B_SplineLU(order, sz, big"0.")
        reslag = fctinter(lag, sz)
        ressplu = fctinter(splu, sz)
        println("$order\t$reslag\t$ressplu")
    end
end
function fctmainsize_gnuplot(szlist, orderlist)

    prec = precision(BigFloat)
    println("# begin gnuplot order=$orderlist julia precision=$prec")
    print("# size")
    for order in orderlist
        print("\tLagrange($order)\tB_Spline($order)")
    end
    println("")
    for sz in szlist
        print("$sz")
        for order in orderlist
            lag = Lagrange(BigFloat, order)
            splu = B_SplineLU(order, sz, big"0.")
            reslag = fctinter(lag, sz)
            ressplu = fctinter(splu, sz)
            print("\t$reslag\t$ressplu")
        end
        println("")
    end
end

setprecision(1024) do
#    fctmain_gnuplot(1001, 51)
    fctmainsize_gnuplot(101:100:10001,[3,11,31])
end

