set title "LipoP predictions for SPy_1798"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1798.eps"
set arrow from 2,5.19094 to 6,5.19094 nohead lt 1 lw 20
set label "SpI" at 7,5.19094
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,5.19094 to 6,5.19094 nohead lt 1 lw 20
set label "SpI" at 7,5.19094
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
26.500000 7.578830
25.500000 4.738290
23.500000 4.236660
22.500000 4.079170
20.500000 3.980803
21.500000 3.577838
24.500000 1.985410
18.500000 1.941800
19.500000 1.546410
e
exit
