set title "LipoP predictions for SPy_0327"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:15]
set y2range [0:18]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0327.eps"
set arrow from 2,13.4022 to 6,13.4022 nohead lt 4 lw 20
set label "TMH" at 7,13.4022
set arrow from 2,4.45676 to 6,4.45676 nohead lt 1 lw 20
set label "SpI" at 7,4.45676
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,13.4022 to 6,13.4022 nohead lt 4 lw 20
set label "TMH" at 7,13.4022
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
42.500000 7.144040
41.500000 3.675369
36.500000 3.312130
43.500000 2.463377
38.500000 1.718830
e
exit
