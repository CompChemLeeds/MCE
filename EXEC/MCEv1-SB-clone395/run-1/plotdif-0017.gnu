set terminal png
set output "PopDiff-0017.png"
set title "Graph of Population Difference"
set ylabel "Population Difference"
plot "normpop-0021.out" u 1:2 t "Norm" w l
