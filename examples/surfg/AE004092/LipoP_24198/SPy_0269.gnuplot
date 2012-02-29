set title "LipoP predictions for SPy_0269"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0269.eps"
set arrow from 2,5.61352 to 6,5.61352 nohead lt 1 lw 20
set label "SpI" at 7,5.61352
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,5.61352 to 6,5.61352 nohead lt 1 lw 20
set label "SpI" at 7,5.61352
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
36.500000 8.552720
31.500000 3.146136
33.500000 1.360460
32.500000 1.276460
30.500000 0.138040
e
exit
