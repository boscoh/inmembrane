set title "LipoP predictions for SPy_0588"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:5]
set y2range [0:8]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0588.eps"
set arrow from 2,2.79408 to 6,2.79408 nohead lt 1 lw 20
set label "SpI" at 7,2.79408
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,2.79408 to 6,2.79408 nohead lt 1 lw 20
set label "SpI" at 7,2.79408
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
23.500000 5.268950
19.500000 3.442007
18.500000 0.882850
22.500000 0.112690
e
exit
