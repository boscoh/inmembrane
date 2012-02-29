set title "LipoP predictions for SPy_0136"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:20]
set y2range [0:23]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0136.eps"
set arrow from 2,16.6463 to 6,16.6463 nohead lt 1 lw 20
set label "SpI" at 7,16.6463
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,16.6463 to 6,16.6463 nohead lt 1 lw 20
set label "SpI" at 7,16.6463
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
27.500000 19.640600
25.500000 10.929340
23.500000 9.774070
22.500000 7.437070
29.500000 7.086440
26.500000 5.244480
20.500000 4.558740
24.500000 4.324780
21.500000 3.560385
e
exit
