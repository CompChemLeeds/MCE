set terminal png
plot "normpop-0047.out" u 1:2 t "Norm" w l
set title "Graph of Norm"
set output "PopDiff-0041.png"
set title "Graph of Population Difference"
plot "normpop-0033.out" u 1:6 t "Real" w l, "" u 1:7 t "Imaginary" w l, "" u 1:8 t "Absolute" w l
