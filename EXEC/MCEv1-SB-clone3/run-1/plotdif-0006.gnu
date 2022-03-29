set terminal png
set output "PopDiff-0006.png"
plot "normpop-0004.out" u 1:3 t "Re(ACF)" w l, "" u 1:4 t "Im(ACF)" w l, "" u 1:5 t "|ACF|" w l
