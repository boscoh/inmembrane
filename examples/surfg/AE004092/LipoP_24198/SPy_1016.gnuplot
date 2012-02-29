set title "LipoP predictions for SPy_1016"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:5]
set y2range [0:8]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1016.eps"
set arrow from 2,3.68845 to 6,3.68845 nohead lt 1 lw 20
set label "SpI" at 7,3.68845
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,3.68845 to 6,3.68845 nohead lt 1 lw 20
set label "SpI" at 7,3.68845
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
23.500000 6.325880
25.500000 3.956614
28.500000 1.150690
24.500000 1.046590
27.500000 0.962880
e
exit
