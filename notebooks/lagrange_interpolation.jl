using Plots

const n = 10
const a = -10.0
const b = 10.0

f( x ) =  1. / (1 .+ x*x)


# +
xp = LinRange(a, b, n)
yp = f.(xp) 

plot(xp, yp)
# -

function lagrange_interpolate( x, xp, yp )
  P = 0.0
  for k in eachindex(yp)
    L = 1.0
    for i in eachindex(xp) 
        L = (i != k) &&  L * (x - xp[i])/(xp[k] - xp[i])
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

xs, ys = interpolation(xp, yp)

plot(xs, ys)


