set title "LipoP predictions for SPy_1225"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:15]
set y2range [0:18]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1225.eps"
set arrow from 2,9.86168 to 6,9.86168 nohead lt 4 lw 20
set label "TMH" at 7,9.86168
set arrow from 2,1.98683 to 6,1.98683 nohead lt 1 lw 20
set label "SpI" at 7,1.98683
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,9.86168 to 6,9.86168 nohead lt 4 lw 20
set label "TMH" at 7,9.86168
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
18.500000 4.650010
27.500000 1.008170
23.500000 0.731090
24.500000 0.308550
e
exit
