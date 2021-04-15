set terminal png font "arial"
set ylabel "Error"
set xlabel "Order"
set title "Error as a function of order for size=1001 and precision=1024 bits"
set logscale y
set grid
set format y "10^{%L}"
# set pointsize 0.6
# set style line 1 lc rgb 'cyan' pt 5   # circle
# set style line 2 lc rgb 'green' pt 5   # triangle
# set yrange [0:100]
# set key left top Left reverse
set output "../out/ord_sz1001_1d.png"
plot [3:51] "ord_sz1001_1d.dat" u 1:2 t "Lagrange" with linespoints pt 1, "" u 1:3 t "Spline" with linespoints pt 1
