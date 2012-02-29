set title "LipoP predictions for SPy_0168"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:25]
set y2range [0:28]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0168.eps"
set arrow from 2,19.7897 to 6,19.7897 nohead lt 1 lw 20
set label "SpI" at 7,19.7897
set arrow from 2,1.32801 to 6,1.32801 nohead lt 4 lw 20
set label "TMH" at 7,1.32801
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,19.7897 to 6,19.7897 nohead lt 1 lw 20
set label "SpI" at 7,19.7897
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
31.500000 22.789700
26.500000 6.339750
33.500000 6.303000
30.500000 4.658690
28.500000 3.850068
25.500000 0.917170
32.500000 0.818130
e
exit
