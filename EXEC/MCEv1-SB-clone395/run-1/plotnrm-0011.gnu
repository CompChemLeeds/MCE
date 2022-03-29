set terminal png
set output "Norm-0011.png"
set title "Graph of Norm"
set xlabel "Time (au)"
set ylabel "Norm"
set ylabel "Extra"
set xlabel "Time (au)"
set ylabel "ACF"
plot "normpop-0007.out" u 1:3 t "Re(ACF)" w l, "" u 1:4 t "Im(ACF)" w l, "" u 1:5 t "|ACF|" w l
