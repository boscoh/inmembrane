set title "LipoP predictions for SPy_0210"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:15]
set y2range [0:18]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0210.eps"
set arrow from 2,12.3566 to 6,12.3566 nohead lt 2 lw 20
set label "SpII" at 7,12.3566
set arrow from 2,1.79029 to 6,1.79029 nohead lt 1 lw 20
set label "SpI" at 7,1.79029
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,12.3566 to 6,12.3566 nohead lt 2 lw 20
set label "SpII" at 7,12.3566
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 2 lw 20, "-" axes x1y2 title "" with impulses lt 1 lw 20
24.500000 15.356600
e
27.500000 3.941949
26.500000 2.607911
29.500000 2.317053
e
exit
