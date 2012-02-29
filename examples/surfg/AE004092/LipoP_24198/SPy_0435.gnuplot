set title "LipoP predictions for SPy_0435"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:15]
set y2range [0:18]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0435.eps"
set arrow from 2,9.83082 to 6,9.83082 nohead lt 1 lw 20
set label "SpI" at 7,9.83082
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,9.83082 to 6,9.83082 nohead lt 1 lw 20
set label "SpI" at 7,9.83082
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
27.500000 12.319660
24.500000 10.910610
22.500000 6.928990
20.500000 6.217800
21.500000 4.817000
26.500000 2.555397
29.500000 2.372027
19.500000 2.248859
15.500000 2.047021
18.500000 0.863550
23.500000 0.800380
e
exit
