set title "LipoP predictions for SPy_1273"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:20]
set y2range [0:23]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1273.eps"
set arrow from 2,15.1377 to 6,15.1377 nohead lt 1 lw 20
set label "SpI" at 7,15.1377
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,15.1377 to 6,15.1377 nohead lt 1 lw 20
set label "SpI" at 7,15.1377
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
28.500000 18.103000
25.500000 12.305400
27.500000 10.703150
23.500000 6.433590
20.500000 3.595700
19.500000 3.569282
30.500000 3.038905
26.500000 1.520210
22.500000 0.633980
e
exit
