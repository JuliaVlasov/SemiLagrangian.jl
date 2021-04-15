set terminal png font "times"
set ylabel "Electric energy"
set xlabel "time"
set title "Electric energy as a function of time with Lagrange order 51"
set logscale y
set grid
set format y "10^{%L}"
# set pointsize 0.6
# set style line 1 lc rgb 'cyan' pt 5   # circle
# set style line 2 lc rgb 'green' pt 5   # triangle
# set yrange [0:100]
# set key left top Left reverse
set output "../out/landau_64_10000_51.png"
plot [0:20] "landau_64_10000_51.txt" u 1:2 t "v=64x64, ε = 0.5" with linespoints pt 1,
 
