set title "LipoP predictions for SPy_0458"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0458.eps"
set arrow from 2,7.01411 to 6,7.01411 nohead lt 1 lw 20
set label "SpI" at 7,7.01411
set arrow from 2,1.84328 to 6,1.84328 nohead lt 4 lw 20
set label "TMH" at 7,1.84328
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,7.01411 to 6,7.01411 nohead lt 1 lw 20
set label "SpI" at 7,7.01411
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
51.500000 9.934920
48.500000 5.057330
47.500000 3.508738
49.500000 2.719536
50.500000 0.308900
e
exit
