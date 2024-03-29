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

# function trace_energy(
#     advd::AdvectionData{T,N,timeopt},
#     t,
# ) where {T,N,timeopt}
#     compute_charge!(advd)
#     compute_elfield!(advd)
#     # clockend(cl_obs,6)
#     # clockbegin(cl_obs,7)
#     elenergy = compute_ee(advd)
#     # clockend(cl_obs,7)
#     # clockbegin(cl_obs,8)
#     kinenergy = compute_ke(advd)
#     # clockend(cl_obs,8)
#     energyall = elenergy + kinenergy
#     return energyall
# end

function landau(advd::AdvectionData, nbdt)

    # global cl_obs
    # clockreset(cl_obs)
    maxdiff = 0
    dt = advd.adv.dt_base
    refel = getenergyall(advd)
    #   @show typeof(refel)
    #    printall(cl_obs)
    #    clockreset(cl_obs)
    maxel = refel
    minel = refel
    for i = 1:nbdt
        while advection!(advd)
        end
        el = getenergyall(advd)
        #       @show typeof(el)
        maxel = max(maxel, el)
        minel = min(minel, el)
        maxdiff = maxel - minel
        # @show i,el,diff, maxdiff
        # flush(stdout)
        # printall(cl_obs)
        # clockreset(cl_obs)
    end
    # printall(cl_obs)
    return maxdiff
end

function landau2(
    t_max::T,
    timeopt,
    dt::T,
    sz,
    interp::AbstractInterpolation,
    type::TimeAlgorithm,
    epsilon::T,
    typeadd,
) where {T}
    nbdt = Int(round(t_max / dt))
    spmin, spmax, nsp = T(0), T(4big(pi)), sz[1]
    vmin, vmax, nv = -T(10), T(10), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    tabst = [([1, 2], 2, 1, false)]

    adv = Advection(
        (mesh_sp, mesh_v),
        [interp, interp],
        dt,
        tabst;
        tab_coef = nosplit(dt),
        timeopt = timeopt,
        timealg = type,
        ordalg = typeadd,
    )

    fct_sp(x) = epsilon * cos(x / 2) + 1
    fct_v(v) = exp(-v^2 / 2) / sqrt(2T(pi))

    lgn_sp = fct_sp.(mesh_sp.points)
    lgn_v = fct_v.(mesh_v.points)

    data = dotprod((lgn_sp, lgn_v))

    resint = sum(data) * (vmax - vmin) / length(data)
    data /= resint

    pvar = getpoissonvar(adv; type = StdPoisson2d, typeadd = typeadd)

    advd = AdvectionData(adv, data, pvar)

    if type == StdABinit
        # tabst_2 = [( [2,1], 1, 1, true, true),( [1,2], 1, 2, true, false) ]
        # adv_2 = Advection(
        # (mesh_sp,mesh_v,), 
        # [interp,interp], 
        # dt,
        # tabst_2,
        # tab_coef=hamsplit_3_11(dt), 
        # timeopt = timeopt)

        # pvar_2 = getpoissonvar(adv_2)

        # advd_2 = AdvectionData(adv_2, data, pvar_2)

        # init_data = map(x->zeros(T,sz), 1:(typeadd-1))
        # for i=1:(typeadd-1)
        #     while advection!(advd_2)
        #     end
        #     pvar.init_data[i] .= advd_2.data
        # end
        for i = 1:(typeadd-1)
            pvar.init_data[i] .= data
        end
    end

    # advdata = Advection1dData(adv, data, pvar)

    return landau(advd, nbdt)
end

function landau1_1(
    t_max::T,
    timeopt,
    dt::T,
    sz,
    interp::AbstractInterpolation,
    tab_coef,
    epsilon::T,
) where {T}
    nbdt = Int(round(t_max / dt))
    spmin, spmax, nsp = T(0), T(4big(pi)), sz[1]
    vmin, vmax, nv = -T(10), T(10), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    tabst = [([2, 1], 1, 1, true), ([1, 2], 1, 2, true)]
    #    tabst = [( [1,2], 1, 1, true),( [2,1], 1, 2, true) ]

    adv = Advection(
        (mesh_sp, mesh_v),
        [interp, interp],
        dt,
        tabst;
        tab_coef = tab_coef,
        timeopt = timeopt,
    )

    fct_sp(x) = epsilon * cos(x / 2) + 1
    fct_v(v) = exp(-v^2 / 2) / sqrt(2T(pi))

    lgn_sp = fct_sp.(mesh_sp.points)
    lgn_v = fct_v.(mesh_v.points)

    data = dotprod((lgn_sp, lgn_v))

    resint = sum(data) * (vmax - vmin) / length(data)
    data /= resint

    pvar = getpoissonvar(adv)

    advd = AdvectionData(adv, data, pvar)

    # advdata = Advection1dData(adv, data, pvar)

    return landau(advd, nbdt)
end

function run_mesure(t_max::T, timeopt, sz, interp, epsilon) where {T}
    # tabsplit = [standardsplit, strangsplit, triplejumpsplit, order6split, hamsplit_3_11]
    # tabsplit = [standardsplit, strangsplit, triplejumpsplit, table2split]
    tabsplit = [standardsplit, strangsplit, hamsplit_3_11]
    # tabtype = [StdPoisson2d, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB, StdAB]
    # tabtype = [NoTimeAlg, ABTimeAlg,ABTimeAlg,ABTimeAlg,ABTimeAlg,ABTimeAlg,ABTimeAlg,]
    # tabtype = [StdRK4,]
    tabtype = []
    tabtypeadd = [0, 2, 3, 4, 5, 6]
    # tabtxtsplit = ["stdsplit", "strangsplit", "triplejumpsplit", "order6split", "fernandosplit"]
    # tabtxt = ["stdsplit", "strangsplit", "fernandosplit", "std2d", "stdAB_2", "stdAB_3", "stdAB_4", "stdAB_5", "stdAB_6",]
    # tabtxt = ["stdsplit", "strangsplit", "fernandosplit", "stdRK4",]
    tabtxt = ["stdsplit", "strangsplit", "fernandosplit"]
    tabnbdt = [
        10,
        20,
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

    res = zeros(Float64, length(tabtxt) + 1, length(tabnbdt))

    res[1, :] .= Float64(t_max) ./ tabnbdt

    lenitc = length(tabsplit)
    for inbdt = 1:length(tabnbdt), itc = 1:length(tabtxt)
        tc = itc <= lenitc ? tabsplit[itc] : missing
        tp = itc > lenitc ? tabtype[itc-lenitc] : missing
        typadd = itc > lenitc ? tabtypeadd[itc-lenitc] : missing
        nbdt = tabnbdt[inbdt]
        dt = t_max / nbdt
        res[itc+1, inbdt] = if itc <= lenitc
            landau1_1(t_max, timeopt, dt, sz, interp, tc(dt), epsilon)
        else
            landau2(t_max, timeopt, dt, sz, interp, tp, epsilon, typadd)
        end
        if timeopt != MPIOpt || MPI.Comm_rank(MPI.COMM_WORLD) == 1
            println("# sz=$sz t_max=$t_max interp=$interp")
            for txt in tabtxt
                print("# $txt\t")
            end
            println("")
            for j = 1:size(res, 2), i = 1:size(res, 1)
                print("$(res[i,j])")
                if i == size(res, 1)
                    print("\n")
                else
                    print("\t")
                end
            end

            println("free memory : $(Sys.free_memory()/2^30)")
            flush(stdout)
        end
    end
end
# to run with mpi here is the command for n process (n must be a number) from SemiLagrangian.jl :
# ./mpiexec.loc n notebooks/run_mesures_poisson.jl
# else run like this :
# julia --project=.  notebooks/run_mesures_poisson.jl

T = Double64
#run_mesure(T(1), NoTimeOpt, (128,128), Lagrange(11,T), T(0.5))
#run_mesure(T(1), MPIOpt, (128,128), Lagrange(11,T), T(0.5))
run_mesure(T(1), SimpleThreadsOpt, (128, 128), Lagrange(11, T), T(0.5))
