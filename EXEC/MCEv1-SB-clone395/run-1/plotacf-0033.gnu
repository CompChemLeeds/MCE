set terminal png
set output "ACF-0033.png"
set title "Graph of Autocorrelation Function"
set output "ACF-0047.png"
set title "Graph of Autocorrelation Function"
set title "Graph of Extra Calculated Quantity"
set xlabel "Time (au)"
set ylabel "ACF"
plot "normpop-0033.out" u 1:3 t "Re(ACF)" w l, "" u 1:4 t "Im(ACF)" w l, "" u 1:5 t "|ACF|" w l
