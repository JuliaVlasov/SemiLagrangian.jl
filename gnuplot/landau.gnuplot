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
plot [0:10] "landau.dat" u 1:2 t "Electric energy" with linespoints pt 1
