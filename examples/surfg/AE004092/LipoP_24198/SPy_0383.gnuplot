set title "LipoP predictions for SPy_0383"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0383.eps"
set arrow from 2,7.06352 to 6,7.06352 nohead lt 4 lw 20
set label "TMH" at 7,7.06352
set arrow from 2,2.59217 to 6,2.59217 nohead lt 1 lw 20
set label "SpI" at 7,2.59217
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,7.06352 to 6,7.06352 nohead lt 4 lw 20
set label "TMH" at 7,7.06352
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
39.500000 4.938830
41.500000 3.155019
37.500000 2.269005
40.500000 0.913930
48.500000 0.077150
e
exit
