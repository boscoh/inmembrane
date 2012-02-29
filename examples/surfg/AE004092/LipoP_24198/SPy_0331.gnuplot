set title "LipoP predictions for SPy_0331"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:20]
set y2range [0:23]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0331.eps"
set arrow from 2,15.4464 to 6,15.4464 nohead lt 4 lw 20
set label "TMH" at 7,15.4464
set arrow from 2,2.02634 to 6,2.02634 nohead lt 1 lw 20
set label "SpI" at 7,2.02634
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,15.4464 to 6,15.4464 nohead lt 4 lw 20
set label "TMH" at 7,15.4464
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
36.500000 3.267670
31.500000 2.868827
32.500000 2.603976
30.500000 2.423426
34.500000 0.556850
e
exit
