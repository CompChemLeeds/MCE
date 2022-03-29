set terminal png
set output "ACF-0004.png"
set title "Graph of Extra Calculated Quantity"
set ylabel "Extra"
set xlabel "Time (au)"
set ylabel "Norm"
plot "normpop-0006.out" u 1:2 t "Norm" w l
set title "Graph of Autocorrelation Function"
set xlabel "Time (au)"
set ylabel "ACF"
plot "normpop-0004.out" u 1:3 t "Re(ACF)" w l, "" u 1:4 t "Im(ACF)" w l, "" u 1:5 t "|ACF|" w l
