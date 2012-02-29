set title "LipoP predictions for SPy_1892"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:5]
set y2range [0:8]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1892.eps"
set arrow from 2,1.07363 to 6,1.07363 nohead lt 4 lw 20
set label "TMH" at 7,1.07363
set arrow from 2,1.04215 to 6,1.04215 nohead lt 1 lw 20
set label "SpI" at 7,1.04215
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,1.07363 to 6,1.07363 nohead lt 4 lw 20
set label "TMH" at 7,1.07363
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
34.500000 2.498190
30.500000 1.345320
32.500000 1.187580
23.500000 0.943380
29.500000 0.506940
31.500000 0.048830
e
exit
