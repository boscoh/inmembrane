set title "LipoP predictions for SPy_0437"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:20]
set y2range [0:23]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0437.eps"
set arrow from 2,16.303 to 6,16.303 nohead lt 1 lw 20
set label "SpI" at 7,16.303
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,16.303 to 6,16.303 nohead lt 1 lw 20
set label "SpI" at 7,16.303
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
22.500000 19.297400
21.500000 10.949160
23.500000 8.060520
24.500000 7.406410
18.500000 5.324530
20.500000 4.628380
19.500000 4.293290
15.500000 1.994320
e
exit
