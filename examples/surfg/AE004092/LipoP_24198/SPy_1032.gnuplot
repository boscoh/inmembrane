set title "LipoP predictions for SPy_1032"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1032.eps"
set arrow from 2,6.31355 to 6,6.31355 nohead lt 1 lw 20
set label "SpI" at 7,6.31355
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,6.31355 to 6,6.31355 nohead lt 1 lw 20
set label "SpI" at 7,6.31355
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
33.500000 8.532050
32.500000 6.881890
30.500000 6.658080
28.500000 4.601800
26.500000 3.882536
29.500000 1.548020
35.500000 0.957980
25.500000 0.451980
e
exit
