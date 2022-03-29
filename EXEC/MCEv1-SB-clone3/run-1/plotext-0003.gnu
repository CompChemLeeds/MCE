set terminal png
set ylabel "Population Difference"
set xlabel "Time (au)"
plot "normpop-0005.out" u 1:13 t "Population Difference" w l, "" u 1:2 t "Norm" w l
