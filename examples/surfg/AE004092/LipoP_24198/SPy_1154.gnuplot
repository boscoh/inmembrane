set title "LipoP predictions for SPy_1154"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1154.eps"
set arrow from 2,5.94269 to 6,5.94269 nohead lt 1 lw 20
set label "SpI" at 7,5.94269
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,5.94269 to 6,5.94269 nohead lt 1 lw 20
set label "SpI" at 7,5.94269
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
34.500000 8.187790
32.500000 6.054470
39.500000 5.692340
38.500000 5.159630
33.500000 4.366270
35.500000 4.238560
40.500000 2.043719
36.500000 0.247670
e
exit
