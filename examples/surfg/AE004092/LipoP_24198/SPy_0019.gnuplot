set title "LipoP predictions for SPy_0019"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:20]
set y2range [0:23]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0019.eps"
set arrow from 2,16.394 to 6,16.394 nohead lt 1 lw 20
set label "SpI" at 7,16.394
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,16.394 to 6,16.394 nohead lt 1 lw 20
set label "SpI" at 7,16.394
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
24.500000 19.268200
23.500000 15.645700
18.500000 11.177830
17.500000 10.298650
19.500000 10.250480
21.500000 9.312900
20.500000 8.828030
22.500000 8.020650
26.500000 5.718430
25.500000 5.349900
27.500000 3.362370
16.500000 1.795730
28.500000 0.161170
e
exit
