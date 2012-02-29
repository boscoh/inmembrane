set title "LipoP predictions for SPy_0230"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0230.eps"
set arrow from 2,5.59408 to 6,5.59408 nohead lt 4 lw 20
set label "TMH" at 7,5.59408
set arrow from 2,5.51575 to 6,5.51575 nohead lt 1 lw 20
set label "SpI" at 7,5.51575
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,5.59408 to 6,5.59408 nohead lt 4 lw 20
set label "TMH" at 7,5.59408
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
33.500000 8.286270
35.500000 4.857530
31.500000 3.955431
39.500000 1.843990
34.500000 1.724500
e
exit
