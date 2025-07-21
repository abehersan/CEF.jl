set terminal wxt enhanced title 'CEF.jl point charges'
set xlabel 'x (Å)'
set ylabel 'y (Å)'
set zlabel 'z (Å)'
set view equal xyz
splot \
    '-' using 1:2:3 with points pointtype 7 pointsize 1.5 lc rgb 'purple' title 'RE ion', \
    '-' using 1:2:3 with points pointtype 7 pointsize 1.5 lc rgb 'orange' title 'PCs', \
    '-' using 1:2:3:4:5:6 with vectors head filled lc rgb 'red' lw 2 title 'a', \
    '-' using 1:2:3:4:5:6 with vectors head filled lc rgb 'green' lw 2 title 'b', \
    '-' using 1:2:3:4:5:6 with vectors head filled lc rgb 'blue' lw 2 title 'c'
0.0 0.0 0.0
e
-1.22069536 1.9039450722745968 1.5667015000000004
2.25928224 0.10518065590449223 1.5667015000000004
-1.0384822399999998 -2.009186142111257 1.5667015000000004
1.2287526400000002 1.908596945051533 -1.5667651999999996
-2.26726976 0.10983252868142844 -1.5667651999999996
1.03851712 -2.0184294737329616 -1.5667651999999996
e
0 0 0 6.976 0.0 0.0
e
0 0 0 -3.488 6.041393216800244 0.0
e
0 0 0 0.0 0.0 19.11
e
pause -1
