set terminal png font "times"
set ylabel "Electric energy"
set xlabel "time"
set title "Electric energy as a function of time with Lagrange interpolation"
set logscale y
set grid
set format y "10^{%L}"
# set pointsize 0.6
# set style line 1 lc rgb 'cyan' pt 5   # circle
# set style line 2 lc rgb 'green' pt 5   # triangle
# set yrange [0:100]
set key left top Left reverse
set output "../out/landau_lagrange.png"
# plot [0:24.5] "landau_64_10000_27.txt" u 1:2 t "sp, v=64x64, ε = 0.5, order=27" with linespoints pt 0, "landau_64_10000_51.txt" u 1:2 t "sp, v=64x64, ε = 0.5, order=51" with linespoints pt 2, "landau_64_10000_101.txt" u 1:2 t "sp, v=64x64, ε = 0.5, order=101" with linespoints pt 3, "landau_128_67.txt" u 1:2 t "sp, v=128x128, ε = 0.5, order=67" with linespoints pt 4,
 
plot [0:80] "trace32_32_5.txt" u 1:2 t "sp=32^2, v=32^2, ε = 0.5, order=5, dt=1/8" with linespoints pt 0, "trace32_128_lag_5.txt" u 1:2 t "sp=32^2, v=128^2, ε = 0.5, order=5, dt=1/8" with linespoints pt 0, "trace32_128_lag_19.txt" u 1:2 t "sp=32^2, v=128^2, ε = 0.5, order=19, dt=1/8" with linespoints pt 0, "trace32_128_lag_19_dt80.txt" u 1:2 t "sp=32^2, v=128^2, ε = 0.5, order=19, dt=1/80" with linespoints pt 0,
 
