set terminal png font "times"
set ylabel "Delta energy"
set xlabel "Delta t"
set title "Delta total energy as a function of delta time"
set logscale x
set logscale y
set grid
set format y "10^{%L}"
# set pointsize 0.6
# set style line 1 lc rgb 'cyan' pt 5   # circle
# set style line 2 lc rgb 'green' pt 5   # triangle
# set yrange [0:100]
set key right top Left reverse
set output "../out/mesures_poisson.png"
# plot [0:24.5] "landau_64_10000_27.txt" u 1:2 t "sp, v=64x64, ε = 0.5, order=27" with linespoints pt 0, "landau_64_10000_51.txt" u 1:2 t "sp, v=64x64, ε = 0.5, order=51" with linespoints pt 2, "landau_64_10000_101.txt" u 1:2 t "sp, v=64x64, ε = 0.5, order=101" with linespoints pt 3, "landau_128_67.txt" u 1:2 t "sp, v=128x128, ε = 0.5, order=67" with linespoints pt 4,
 
plot [0.5:0.001] "mesures_poisson.txt" u 1:2 t "no split" with linespoints pt 1,  "mesures_poisson.txt" u 1:3 t "Strang split" with linespoints pt 2,  "mesures_poisson.txt" u 1:4 t "triple jump split" with linespoints pt 3,  "mesures_poisson.txt" u 1:5 t "order6 split" with linespoints pt 4 
 
