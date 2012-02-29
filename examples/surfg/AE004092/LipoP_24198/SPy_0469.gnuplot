set title "LipoP predictions for SPy_0469"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:20]
set y2range [0:23]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0469.eps"
set arrow from 2,14.6664 to 6,14.6664 nohead lt 1 lw 20
set label "SpI" at 7,14.6664
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,14.6664 to 6,14.6664 nohead lt 1 lw 20
set label "SpI" at 7,14.6664
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
25.500000 17.442000
23.500000 14.761300
27.500000 9.605740
21.500000 9.567700
24.500000 9.038170
26.500000 6.165950
22.500000 5.890740
28.500000 1.508700
e
exit
