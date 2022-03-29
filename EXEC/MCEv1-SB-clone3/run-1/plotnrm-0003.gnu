set terminal png
set output "Norm-0003.png"
set title "Graph of Norm"
set xlabel "Time (au)"
set title "Graph of Population Difference"
set ylabel "Norm"
plot "normpop-0003.out" u 1:2 t "Norm" w l
