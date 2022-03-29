set terminal png
set output "Extra-0059.png"
set title "Graph of Extra Calculated Quantity"
set ylabel "Extra"
set xlabel "Time (au)"
plot "normpop-0059.out" u 1:6 t "Real" w l, "" u 1:7 t "Imaginary" w l, "" u 1:8 t "Absolute" w l
