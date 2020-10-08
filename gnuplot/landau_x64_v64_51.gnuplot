set terminal png font "times"
set ylabel "Electric energy"
set xlabel "time"
set title "Electric energy as a function of time with Lagrange order 41"
set logscale y
set grid
set format y "10^{%L}"
# set pointsize 0.6
# set style line 1 lc rgb 'cyan' pt 5   # circle
# set style line 2 lc rgb 'green' pt 5   # triangle
# set yrange [0:100]
# set key left top Left reverse
set output "../out/landau_x64_v64_51.png"
plot [0:200] "landau_x64_v64_51_200.dat" u 1:2 t "v=64, ε = 10^{-3}" with linespoints pt 1,
 
