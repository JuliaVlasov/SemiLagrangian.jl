# -*- coding: utf-8 -*-
using Plots

n = 21
a = -5.0
b = 5.0

f( x ) =  exp( - x^2 / 2)


# +
#xp = [0.5 * (a+b) + 0.5 * (b-a) * cos((2k-1)π/ 2n) for k in 1:n]
xp = LinRange(a, b, n)
yp = f.(xp) 

plot(xp, yp)
# -

function lagrange_interpolate( x, xp, yp )
    P = 0.0
    for k in eachindex(yp)
        L = 1.0
        for i in eachindex(xp) 
            if (i != k)  L = L * (x - xp[i])/(xp[k] - xp[i]) end
        end
        P = P + L * yp[k]
    end
    P
end

function interpolation(xp, yp)
 
   
    println(" Nombre de points d'appui : $n")
    println(" Borne inferieure : $a")
    println(" Borne superieure : $b")


    println("\n data1 Interpolation de Lagrange ");
    println("\n points regulierement espaces et racines de tchebishev ");
    dx = (b-a)/n;
    
    
    ns = 100
    xs = LinRange(a, b, ns)
    ys = zeros(Float64, ns)
    
    for i in eachindex(ys)
    
      ys[i] = lagrange_interpolate(xs[i], xp, yp)

    end
    
    xs, ys
end

xs, ys = interpolation(xp, yp);

plot(xs, ys)
scatter!( xp, yp)


