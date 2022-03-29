set terminal png
set output "Extra-0021.png"
set title "Graph of Extra Calculated Quantity"
set title "Graph of Extra Calculated Quantity"
set ylabel "Extra"
plot "normpop-0029.out" u 1:13 t "Population Difference" w l, "" u 1:2 t "Norm" w l
set ylabel "Extra"
set xlabel "Time (au)"
plot "normpop-0021.out" u 1:6 t "Real" w l, "" u 1:7 t "Imaginary" w l, "" u 1:8 t "Absolute" w l
