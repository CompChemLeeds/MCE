set terminal png
set output "ACF-0011.png"
set title "Graph of Autocorrelation Function"
set xlabel "Time (au)"
set ylabel "ACF"
plot "normpop-0011.out" u 1:3 t "Re(ACF)" w l, "" u 1:4 t "Im(ACF)" w l, "" u 1:5 t "|ACF|" w l
