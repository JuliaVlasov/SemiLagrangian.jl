function aff_graph(ind, fp)
    println("ind=$ind begin")
    for (i, val) in enumerate(fp)
        println("$i\t$val")
    end
    println("ind=$ind end")
end
function test_dirac(bsp::InterpolationType, order, len, nb, modval )
    lag = Lagrange(BigFloat, order; iscirc=true)
#        println("lag=$lag")
    n = len
    # whennot circuar coef != 1 set the function unperiodic
    fp = zeros(BigFloat,n)
    fp[div(n,2)] = big"1.0"
    fi = zeros(BigFloat, n)
    value = big"0.30529810681113334445566767713091"
    normax=0.0
    aff_graph(1,fp)
    for i=1:nb
        fi .= fp
  interpolate!(missing,fp,fi, value, bsp)
        if i%modval == 0
            aff_graph(i,fp)
        end
    end
end
test_dirac(31,100, 500, 50)
