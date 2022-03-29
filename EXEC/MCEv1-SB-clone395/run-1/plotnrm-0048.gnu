set terminal png
set output "Norm-0048.png"
set title "Graph of Norm"
set xlabel "Time (au)"
set ylabel "Norm"
plot "normpop-0048.out" u 1:2 t "Norm" w l
