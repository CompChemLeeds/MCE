set terminal png
set output "Extra-0033.png"
set title "Graph of Extra Calculated Quantity"
set ylabel "Extra"
set xlabel "Time (au)"
plot "normpop-0035.out" u 1:3 t "Re(ACF)" w l, "" u 1:4 t "Im(ACF)" w l, "" u 1:5 t "|ACF|" w l
