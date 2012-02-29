set title "LipoP predictions for SPy_2009"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_2009.eps"
set arrow from 2,8.77151 to 6,8.77151 nohead lt 1 lw 20
set label "SpI" at 7,8.77151
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,-1.62429 to 6,-1.62429 nohead lt 4 lw 20
set label "TMH" at 7,-1.62429
set arrow from 2,8.77151 to 6,8.77151 nohead lt 1 lw 20
set label "SpI" at 7,8.77151
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
37.500000 11.759430
36.500000 3.345106
32.500000 3.338478
35.500000 2.647460
33.500000 0.167010
e
exit
