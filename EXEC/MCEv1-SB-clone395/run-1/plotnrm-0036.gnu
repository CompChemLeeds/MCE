set terminal png
set output "Norm-0036.png"
set title "Graph of Norm"
set xlabel "Time (au)"
set ylabel "Norm"
plot "normpop-0036.out" u 1:2 t "Norm" w l
