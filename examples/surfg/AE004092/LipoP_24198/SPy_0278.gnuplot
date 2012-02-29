set title "LipoP predictions for SPy_0278"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0278.eps"
set arrow from 2,4.87102 to 6,4.87102 nohead lt 1 lw 20
set label "SpI" at 7,4.87102
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,4.87102 to 6,4.87102 nohead lt 1 lw 20
set label "SpI" at 7,4.87102
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
21.500000 7.315370
22.500000 6.029210
23.500000 2.712600
20.500000 1.131300
e
exit
