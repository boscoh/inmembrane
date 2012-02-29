set title "LipoP predictions for SPy_0836"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0836.eps"
set arrow from 2,7.80746 to 6,7.80746 nohead lt 1 lw 20
set label "SpI" at 7,7.80746
set arrow from 2,3.60734 to 6,3.60734 nohead lt 4 lw 20
set label "TMH" at 7,3.60734
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,7.80746 to 6,7.80746 nohead lt 1 lw 20
set label "SpI" at 7,7.80746
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
33.500000 10.440080
34.500000 8.601250
31.500000 2.616083
29.500000 1.949800
32.500000 1.133950
e
exit
