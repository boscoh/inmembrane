set title "LipoP predictions for SPy_1105"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:5]
set y2range [0:8]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1105.eps"
set arrow from 2,3.96621 to 6,3.96621 nohead lt 1 lw 20
set label "SpI" at 7,3.96621
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,-2.77435 to 6,-2.77435 nohead lt 4 lw 20
set label "TMH" at 7,-2.77435
set arrow from 2,3.96621 to 6,3.96621 nohead lt 1 lw 20
set label "SpI" at 7,3.96621
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
27.500000 6.016240
24.500000 4.949590
22.500000 4.505600
25.500000 2.470421
e
exit
