using LinearAlgebra
using DoubleFloats
# using SemiLagMPI
# using MPI
using SemiLagrangian

function printout(
    advd::AdvectionData{T,N,timeopt},
    str,
) where {T,N,timeopt}
    if timeopt != MPIOpt || advd.adv.mpid.ind == 1
        println(str)
    end
end
printout(str) = println(str)

function trace_energy(
    advd::AdvectionData{T,N,timeopt},
    t,
) where {T,N,timeopt}
    compute_charge!(advd)
    compute_elfield!(advd)
    # clockend(cl_obs,6)
    # clockbegin(cl_obs,7)
    elenergy = compute_ee(advd)
    # clockend(cl_obs,7)
    # clockbegin(cl_obs,8)
    kinenergy = compute_ke(advd)
    # clockend(cl_obs,8)
    energyall = elenergy + kinenergy
    return energyall
end

function landau(advd::AdvectionData, nbdt)

    # global cl_obs
    # clockreset(cl_obs)
    maxdiff = 0
    dt = advd.adv.dt_base
    refel = trace_energy(advd, 0*dt)
 #   @show typeof(refel)
    #    printall(cl_obs)
    #    clockreset(cl_obs)
    for i = 1:nbdt
        while advection!(advd)
        end
        el = trace_energy(advd, i * dt)
 #       @show typeof(el)
        maxdiff = max(maxdiff,abs(refel - el))
        # printall(cl_obs)
        # clockreset(cl_obs)
    end
    # printall(cl_obs)
    return maxdiff
end
function landau1_1(
    t_max::T,
    timeopt,
    dt::T,
    sz,
    interp::AbstractInterpolation,
    tab_coef,
    epsilon::T
) where {T}
    nbdt = Int(round(t_max/dt))
    spmin, spmax, nsp = T(0), T(4big(pi)), sz[1]
    vmin, vmax, nv = -T(10), T(10), sz[2]

    mesh_sp = UniformMesh(spmin, spmax, nsp)
    mesh_v = UniformMesh(vmin, vmax, nv)

    tabst = [( [2,1], 1, 1, true, true),( [1,2], 1, 2, true, false) ]

    adv = Advection(
        (mesh_sp,mesh_v,), 
        [interp,interp], 
        dt,
        tabst,
        tab_coef=tab_coef, 
        timeopt = timeopt)

    fct_sp(x) = epsilon * cos(x / 2) + 1
    fct_v(v) = exp(-v^2 / 2) / sqrt(2T(pi))

    lgn_sp = fct_sp.(mesh_sp.points)
    lgn_v = fct_v.(mesh_v.points)

    data = dotprod((lgn_sp, lgn_v))

    resint = sum(data)*(vmax-vmin)/length(data)
    data /= resint
    
    pvar = getpoissonvar(adv)

    advd = AdvectionData(adv, data, pvar)

    # advdata = Advection1dData(adv, data, pvar)

    return landau(advd, nbdt)
end

function run_mesure(    
    t_max::T,
    timeopt,
    sz,
    interp,
    epsilon
) where{T}
# tabsplit = [standardsplit, strangsplit, triplejumpsplit, order6split, hamsplit_3_11]
tabsplit = [standardsplit, strangsplit, triplejumpsplit, hamsplit_3_11, ymsplit]
# tabtxtsplit = ["stdsplit", "strangsplit", "triplejumpsplit", "order6split", "fernandosplit"]
tabtxtsplit = ["stdsplit", "strangsplit", "triplejumpsplit", "fernandosplit", "ymsplit"]
tabnbdt = [10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000]

    res = zeros(Float64, length(tabsplit)+1, length(tabnbdt))

    res[1,:] .= Float64(t_max) ./ tabnbdt 

#        if MPI.Comm_rank(MPI.COMM_WORLD) == 1
#        end


    for inbdt=1:length(tabnbdt), itc=1:length(tabsplit)
        tc = tabsplit[itc]
        nbdt = tabnbdt[inbdt]
        dt = t_max/nbdt
        res[itc+1, inbdt] = landau1_1(t_max, timeopt, dt, sz, interp, tc(dt), epsilon)
#        if MPI.Comm_rank(MPI.COMM_WORLD) == 1
            println("sz=$sz t_max=$t_max interp=$interp")
            for txt in tabtxtsplit
                print("$txt\t")
            end
            println("")
            for j=1:size(res,2),i=1:size(res,1)
                        print("$(res[i,j])")
                if i == size(res,1)
                    print("\n")
                else
                    print("\t")
                end
            end
            GC.gc()
	    println("free memory : $(Sys.free_memory()/2^30)")
#        end
    end
end
T=Double64
run_mesure(T(1), NoTimeOpt, (256,256), Lagrange(19,T), T(0.5))




