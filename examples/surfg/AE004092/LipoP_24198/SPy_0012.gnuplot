set title "LipoP predictions for SPy_0012"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:15]
set y2range [0:18]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0012.eps"
set arrow from 2,11.523 to 6,11.523 nohead lt 1 lw 20
set label "SpI" at 7,11.523
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,11.523 to 6,11.523 nohead lt 1 lw 20
set label "SpI" at 7,11.523
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
21.500000 14.519900
18.500000 3.910420
20.500000 3.746533
23.500000 3.563627
16.500000 1.762560
19.500000 1.294420
22.500000 1.125290
24.500000 1.076080
e
exit
