set terminal png
set output "PopDiff-0012.png"
set title "Graph of Population Difference"
set ylabel "Population Difference"
set xlabel "Time (au)"
set xlabel "Time (au)"
set output "Norm-0010.png"
set output "Extra-0007.png"
set title "Graph of Norm"
set title "Graph of Extra Calculated Quantity"
set xlabel "Time (au)"
set ylabel "Extra"
plot "normpop-0009.out" u 1:13 t "Population Difference" w l, "" u 1:2 t "Norm" w l
plot "normpop-0011.out" u 1:6 t "Real" w l, "" u 1:7 t "Imaginary" w l, "" u 1:8 t "Absolute" w l
set ylabel "Norm"
