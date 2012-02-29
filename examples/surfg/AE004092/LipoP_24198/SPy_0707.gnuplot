set title "LipoP predictions for SPy_0707"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:5]
set y2range [0:8]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0707.eps"
set arrow from 2,4.0409 to 6,4.0409 nohead lt 1 lw 20
set label "SpI" at 7,4.0409
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,-1.67632 to 6,-1.67632 nohead lt 4 lw 20
set label "TMH" at 7,-1.67632
set arrow from 2,4.0409 to 6,4.0409 nohead lt 1 lw 20
set label "SpI" at 7,4.0409
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
25.500000 6.756080
27.500000 4.060210
30.500000 1.818210
28.500000 0.768500
e
exit
