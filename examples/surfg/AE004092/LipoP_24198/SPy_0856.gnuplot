set title "LipoP predictions for SPy_0856"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:5]
set y2range [0:8]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0856.eps"
set arrow from 2,1.40658 to 6,1.40658 nohead lt 1 lw 20
set label "SpI" at 7,1.40658
set arrow from 2,-0.00988251 to 6,-0.00988251 nohead lt 4 lw 20
set label "TMH" at 7,-0.00988251
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,1.40658 to 6,1.40658 nohead lt 1 lw 20
set label "SpI" at 7,1.40658
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
36.500000 3.648228
30.500000 1.959640
35.500000 0.730230
34.500000 0.581270
e
exit