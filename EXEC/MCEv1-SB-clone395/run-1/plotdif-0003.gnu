set terminal png
set output "PopDiff-0003.png"
set title "Graph of Population Difference"
set ylabel "Population Difference"
set output "Extra-0010.png"
plot "normpop-0011.out" u 1:13 t "Population Difference" w l, "" u 1:2 t "Norm" w l
