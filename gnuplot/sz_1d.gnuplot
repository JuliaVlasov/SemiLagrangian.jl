set terminal png font "arial"
set ylabel "Error"
set xlabel "Size"
set title "Error as a function of size with precision=1024 bits"
set logscale y
set format y "10^{%L}"
# set pointsize 0.6
# set style line 1 lc rgb 'cyan' pt 5   # circle
# set style line 2 lc rgb 'green' pt 5   # triangle
# set yrange [0:100]
# set key right bottom Left reverse
# set key right bottom
set key at graph 1.0, 0.6 # Left reverse
set output "../out/sz_1d.png"
plot [101:10001] "sz_1d.dat" u 1:2 t "Lagrange(order=3)" with linespoints pt 1, "" u 1:3 t "Spline(order=3)" with linespoints pt 1, "" u 1:4 t "Lagrange(order=11)" with linespoints pt 1, "" u 1:5 t "Spline(order=11)" with linespoints pt 1, "" u 1:6 t "Lagrange(order=31)" with linespoints pt 1, "" u 1:7 t "Spline(order=31)" with linespoints pt 1
