set terminal png
set output "Norm-0007.png"
set title "Graph of Norm"
plot "normpop-0012.out" u 1:6 t "Real" w l, "" u 1:7 t "Imaginary" w l, "" u 1:8 t "Absolute" w l
plot "normpop-0011.out" u 1:2 t "Norm" w l
