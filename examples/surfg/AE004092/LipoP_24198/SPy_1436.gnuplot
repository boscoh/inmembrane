set title "LipoP predictions for SPy_1436"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1436.eps"
set arrow from 2,7.95301 to 6,7.95301 nohead lt 1 lw 20
set label "SpI" at 7,7.95301
set arrow from 2,5.39582 to 6,5.39582 nohead lt 4 lw 20
set label "TMH" at 7,5.39582
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,7.95301 to 6,7.95301 nohead lt 1 lw 20
set label "SpI" at 7,7.95301
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
32.500000 10.884090
31.500000 5.836840
33.500000 4.341810
29.500000 2.501666
28.500000 2.023339
30.500000 1.774810
34.500000 0.491120
e
exit
