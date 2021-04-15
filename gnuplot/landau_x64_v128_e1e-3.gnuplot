set terminal png font "arial"
set ylabel "Electric energy"
set xlabel "time"
set title "Error as a function of order for size=1001 and precision=1024 bits"
set logscale y
set grid
set format y "10^{%L}"
# set pointsize 0.6
# set style line 1 lc rgb 'cyan' pt 5   # circle
# set style line 2 lc rgb 'green' pt 5   # triangle
# set yrange [0:100]
# set key left top Left reverse
set output "../out/landau_x64_v128_e1e-3.png
plot [0:100] "landau_x64_v128_e1e-3.dat" u 1:2 with linespoints pt 1
