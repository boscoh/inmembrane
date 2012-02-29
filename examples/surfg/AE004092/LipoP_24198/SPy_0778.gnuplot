set title "LipoP predictions for SPy_0778"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:15]
set y2range [0:18]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0778.eps"
set arrow from 2,13.7306 to 6,13.7306 nohead lt 2 lw 20
set label "SpII" at 7,13.7306
set arrow from 2,3.20182 to 6,3.20182 nohead lt 1 lw 20
set label "SpI" at 7,3.20182
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,13.7306 to 6,13.7306 nohead lt 2 lw 20
set label "SpII" at 7,13.7306
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 2 lw 20, "-" axes x1y2 title "" with impulses lt 1 lw 20
21.500000 16.730600
e
23.500000 6.047500
24.500000 1.618780
25.500000 1.463180
e
exit
