set terminal png
set output "PopDiff-0025.png"
set title "Graph of Population Difference"
set ylabel "Population Difference"
set xlabel "Time (au)"
plot "normpop-0025.out" u 1:13 t "Population Difference" w l, "" u 1:2 t "Norm" w l
set xlabel "Time (au)"
plot "normpop-0013.out" u 1:6 t "Real" w l, "" u 1:7 t "Imaginary" w l, "" u 1:8 t "Absolute" w l
