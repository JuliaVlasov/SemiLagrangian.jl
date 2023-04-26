# # Dirac interpolation

using SemiLagrangian

function aff_graph(ind, fp)
    println("ind=$ind begin")
    minval = 0
    maxval = 0
    for (i, val) in enumerate(fp)
        #        println("$i\t$val")
        minval = min(val, minval)
        maxval = max(val, maxval)
    end
    println("minval=$minval maxval=$maxval")
    println("ind=$ind end")
end

function test_dirac(bsp, order, len, nb, modval)
    #        println("lag=$lag")
    n = len
    # whennot circuar coef != 1 set the function unperiodic
    fp = zeros(BigFloat, n)
    fp[div(n, 2)] = big"1.0"
    fi = zeros(BigFloat, n)
    value = big"0.30529810681113334445566767713091"
    normax = 0.0
    aff_graph(1, fp)
    for i = 1:nb
        fi .= fp
        interpolate!(missing, fp, fi, value, bsp)
        if i % modval == 0
            aff_graph(i, fp)
        end
    end
end

?Lagrange

function test_dirac_lagrange(order, len, nb, modval)
    lag = Lagrange(BigFloat, CircEdge, order)
    return test_dirac(lag, order, len, nb, modval)
end
test_dirac_lagrange(31, 256, 500, 50)


function test_dirac_splu(order, len, nb, modval)
    bsp = B_SplineLU(order, len, BigFloat; iscirc = true)
    return test_dirac(bsp, order, len, nb, modval)
end

function test_dirac_spfft(order, len, nb, modval)
    bsp = B_SplineFFT(order, len, BigFloat)
    return test_dirac(bsp, order, len, nb, modval)
end

test_dirac_splu(31, 256, 500, 50)
test_dirac_spfft(31, 256, 500, 50)


