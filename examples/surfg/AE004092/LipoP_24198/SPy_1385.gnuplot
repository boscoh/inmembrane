set title "LipoP predictions for SPy_1385"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:15]
set y2range [0:18]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1385.eps"
set arrow from 2,9.63533 to 6,9.63533 nohead lt 4 lw 20
set label "TMH" at 7,9.63533
set arrow from 2,4.34001 to 6,4.34001 nohead lt 1 lw 20
set label "SpI" at 7,4.34001
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,9.63533 to 6,9.63533 nohead lt 4 lw 20
set label "TMH" at 7,9.63533
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
31.500000 6.813870
33.500000 5.109100
34.500000 3.256415
36.500000 0.872660
29.500000 0.458680
e
exit
