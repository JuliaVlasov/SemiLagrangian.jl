using LinearAlgebra
using DoubleFloats
# using SemiLagMPI
using MPI
using SemiLagrangian

function printout(advd::AdvectionData{T,N,timeopt}, str) where {T,N,timeopt}
    if timeopt != MPIOpt || advd.adv.mpid.ind == 1
        println(str)
    end
end
printout(str) = println(str)

function pois2(
    t_max::T,
    timeopt,
    dt::T,
    sz,
    interp::AbstractInterpolation,
    typealg::TimeAlgorithm,
    ordalg::Int,
    epsilon::T,
) where {T}
    nbdt = Int(round(t_max / dt))
    spmin, spmax, nsp = T(0), T(4big(pi)), sz[1]
    vmin, vmax, nv = T(-10), T(10), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    tabst = [([1, 2], 2, 1, false, false)]

    adv = Advection(
        (mesh_sp, mesh_v),
        [interp, interp],
        dt,
        tabst;
        tab_coef = nosplit(dt),
        timeopt = timeopt,
        timealg = typealg,
        ordalg = ordalg,
    )

    fct_sp(x) = epsilon * cos(x / 2) + 1
    fct_v(v) = exp(-v^2 / 2) / sqrt(2T(pi))

    lgn_sp = fct_sp.(mesh_sp.points)
    lgn_v = fct_v.(mesh_v.points)

    data = dotprod((lgn_sp, lgn_v))

    resint = sum(data) * (vmax - vmin) / length(data)
    data /= resint

    pvar = getpoissonvar(adv; type = StdPoisson2d)

    result = Array{T,2}[]
    initdatas = missing
    if typealg == ABTimeAlg_init2
        nbfordt = Int(floor(sqrt(nbdt))) + ordalg + 2
        newnb = nbfordt * (3ordalg - 1)
        newdt = dt / nbfordt
        newt_max = dt * (3ordalg - 1)
        newtypealg = ordalg == 2 ? NoTimeAlg : typealg
        newordalg = ordalg == 2 ? 0 : (ordalg - 1)
        res, _ = pois2(newt_max, timeopt, newdt, sz, interp, newtypealg, newordalg, epsilon)
        initdatas = res[nbfordt:nbfordt:end]
        append!(result, initdatas)
    end

    advd = AdvectionData(adv, data, pvar; initdatas = initdatas)

    borne_t = t_max - dt / 2

    maxelall = minelall = getenergyall(advd)

    while advd.time_cur < borne_t
        while advection!(advd)
        end
        push!(result, copy(advd.data))
        el = getenergyall(advd)
        maxelall = max(el, maxelall)
        minelall = min(el, minelall)
    end
    return result, (maxelall - minelall)
end

function run_mesure(t_max::T, timeopt, sz, interp, epsilon::T) where {T}
    tabtypealg = [
        NoTimeAlg,
        ABTimeAlg_ip,
        ABTimeAlg_init2,
        ABTimeAlg_init2,
        ABTimeAlg_init2,
        ABTimeAlg_init2,
    ]
    # ta

    tabordalg = [0, 2, 2, 3, 4, 5]
    # tabtxtsplit = ["stdsplit", "strangsplit", "triplejumpsplit", "order6split", "fernandosplit"]
    # tabtxt = ["stdsplit", "strangsplit", "fernandosplit", "std2d", "stdAB_2", "stdAB_3", "stdAB_4", "stdAB_5", "stdAB_6",]
    # tabtxt = ["stdsplit", "strangsplit", "fernandosplit", "stdRK4",]
    tabtxt = ["NO", "AB_ip2", "AB_init2", "AB_init3", "AB_init4", "AB_init5"]
    tabnbdt = [
        50,
        100,
        200,
        500,
        1000,
        2000,
        5000,
        10000,
        20000,
        50000,
        100000,
        200000,
        500000,
        1000000,
    ]

    #        if MPI.Comm_rank(MPI.COMM_WORLD) == 1
    #        end

    resdata = [zeros(T, sz) for i in CartesianIndices((length(tabnbdt), length(tabtxt)))]
    resen = [zero(T) for i in CartesianIndices((length(tabnbdt), length(tabtxt)))]

    lastind = zeros(Int, length(tabtxt))
    t_loc = time_ns()
    for inbdt = 1:length(tabnbdt), itc = 1:length(tabtxt)
        lastind[itc] = inbdt
        tabres, en = pois2(
            t_max,
            timeopt,
            t_max / tabnbdt[inbdt],
            sz,
            interp,
            tabtypealg[itc],
            tabordalg[itc],
            epsilon,
        )
        resdata[inbdt, itc] .= tabres[end]
        resen[inbdt, itc] = en
        if MPI.Comm_rank(MPI.COMM_WORLD) == 1
            nb = MPI.Comm_size(MPI.COMM_WORLD)
            t_now = time_ns()
            t = (t_now - t_loc) * 1e-9
            t_loc = t_now
            for k = 1:length(tabtxt)
                if lastind[k] > 0 && (k != 1 || lastind[k] > 1)
                    println(
                        "# sz=$sz t_max=$t_max interp=$interp ref=$(tabtxt[k]) nb=$nb t=$t\n# t",
                    )
                    for txt in tabtxt
                        print("\t$txt")
                    end
                    println("")
                    for j = 1:length(tabnbdt), i = 0:length(tabtxt)
                        if i == 0
                            res = Float32(t_max / tabnbdt[j])
                            print("$res")
                        else
                            res = Float64(
                                if (j < lastind[k] || (j == lastind[k] && i < k))
                                    norm(resdata[j, i] - resdata[lastind[k], k])
                                else
                                    0
                                end,
                            )
                            en = resen[j, i]
                            print("$res($en)")
                        end
                        if i == length(tabtxt)
                            print("\n")
                        else
                            print("\t")
                        end
                    end
                end
            end

            println("free memory : $(Sys.free_memory()/2^30)")
            flush(stdout)
        end
    end
end
T = Double64
run_mesure(T(1), MPIOpt, (128, 128), Lagrange(11, T), T(0.5))
