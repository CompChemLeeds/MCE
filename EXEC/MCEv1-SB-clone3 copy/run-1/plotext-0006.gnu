set terminal png
set title "Graph of Population Difference"
set ylabel "Population Difference"
plot "normpop-0003.out" u 1:3 t "Re(ACF)" w l, "" u 1:4 t "Im(ACF)" w l, "" u 1:5 t "|ACF|" w l
