set terminal png
plot "normpop-0012.out" u 1:2 t "Norm" w l
set output "Extra-0009.png"
set title "Graph of Extra Calculated Quantity"
set title "Graph of Autocorrelation Function"
set xlabel "Time (au)"
set ylabel "ACF"
