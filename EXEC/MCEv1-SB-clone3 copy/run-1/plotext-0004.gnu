set terminal png
set xlabel "Time (au)"
plot "normpop-0006.out" u 1:6 t "Real" w l, "" u 1:7 t "Imaginary" w l, "" u 1:8 t "Absolute" w l
