set terminal png
set output "PopDiff-0005.png"
set title "Graph of Population Difference"
set ylabel "Population Difference"
set xlabel "Time (au)"
set ylabel "Norm"
plot "normpop-0004.out" u 1:2 t "Norm" w l
plot "normpop-0006.out" u 1:13 t "Population Difference" w l, "" u 1:2 t "Norm" w l
