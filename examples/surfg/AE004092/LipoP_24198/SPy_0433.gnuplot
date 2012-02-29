set title "LipoP predictions for SPy_0433"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:25]
set y2range [0:28]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0433.eps"
set arrow from 2,19.7965 to 6,19.7965 nohead lt 1 lw 20
set label "SpI" at 7,19.7965
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,19.7965 to 6,19.7965 nohead lt 1 lw 20
set label "SpI" at 7,19.7965
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
22.500000 22.790900
21.500000 14.442600
23.500000 11.553990
24.500000 10.899870
18.500000 8.818000
20.500000 8.121850
19.500000 7.786760
15.500000 5.487790
25.500000 3.299647
17.500000 1.327300
e
exit
