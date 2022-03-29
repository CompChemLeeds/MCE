set terminal png
set ylabel "Extra"
set xlabel "Time (au)"
plot "normpop-0041.out" u 1:6 t "Real" w l, "" u 1:7 t "Imaginary" w l, "" u 1:8 t "Absolute" w l
set output "Norm-0047.png"
set title "Graph of Norm"
