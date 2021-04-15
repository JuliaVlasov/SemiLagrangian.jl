set terminal png font "arial"
set ylabel "Electric energy"
set xlabel "Time"
set title "Electric energy as a function of time"
set logscale y
set grid
set format y "10^{%L}"
# set pointsize 0.6
# set style line 1 lc rgb 'cyan' pt 5   # circle
# set style line 2 lc rgb 'green' pt 5   # triangle
# set yrange [0:100]
# set key left top Left reverse
set output "../out/landau.png"
plot [0:100] "landau1d_128_21.txt" u 1:2 t "x=128, v=128 dt=0.01 lagrange 21, Float64" with linespoints pt 0, "landau_128_21_Double.txt" u 1:2 t "x=128, v=128 dt=0.01 lagrange 21 Float128" with linespoints pt 0, "landau_128_5.txt" u 1:2 t "x=128, v=128 dt=0.01 lagrange 5 Float64" with linespoints pt 0, "landau_256_19_Double.txt" u 1:2 t "x=256, v=256 dt=0.01 lagrange 19 Float128" with linespoints pt 0,
