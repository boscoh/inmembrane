set title "LipoP predictions for SPy_0797"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:15]
set y2range [0:18]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0797.eps"
set arrow from 2,13.0692 to 6,13.0692 nohead lt 4 lw 20
set label "TMH" at 7,13.0692
set arrow from 2,2.58198 to 6,2.58198 nohead lt 1 lw 20
set label "SpI" at 7,2.58198
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,13.0692 to 6,13.0692 nohead lt 4 lw 20
set label "TMH" at 7,13.0692
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
39.500000 4.899670
40.500000 3.601450
38.500000 0.435810
36.500000 0.324820
e
exit
