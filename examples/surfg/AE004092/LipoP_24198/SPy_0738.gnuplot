set title "LipoP predictions for SPy_0738"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:5]
set y2range [0:8]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0738.eps"
set arrow from 2,3.63572 to 6,3.63572 nohead lt 1 lw 20
set label "SpI" at 7,3.63572
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,3.63572 to 6,3.63572 nohead lt 1 lw 20
set label "SpI" at 7,3.63572
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
19.500000 6.381370
21.500000 2.192880
36.500000 2.049781
38.500000 0.987310
16.500000 0.595140
14.500000 0.225470
e
exit
