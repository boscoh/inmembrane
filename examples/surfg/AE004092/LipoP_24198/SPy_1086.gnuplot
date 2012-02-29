set title "LipoP predictions for SPy_1086"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1086.eps"
set arrow from 2,6.42548 to 6,6.42548 nohead lt 4 lw 20
set label "TMH" at 7,6.42548
set arrow from 2,4.26337 to 6,4.26337 nohead lt 1 lw 20
set label "SpI" at 7,4.26337
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,6.42548 to 6,6.42548 nohead lt 4 lw 20
set label "TMH" at 7,6.42548
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
39.500000 7.128400
37.500000 2.964294
40.500000 1.803350
36.500000 0.703470
e
exit
