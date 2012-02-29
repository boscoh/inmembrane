set title "LipoP predictions for SPy_0128"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:20]
set y2range [0:23]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0128.eps"
set arrow from 2,16.7355 to 6,16.7355 nohead lt 1 lw 20
set label "SpI" at 7,16.7355
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,16.7355 to 6,16.7355 nohead lt 1 lw 20
set label "SpI" at 7,16.7355
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
23.500000 19.696100
17.500000 14.248100
20.500000 10.679680
22.500000 10.571940
21.500000 9.012380
19.500000 8.179670
18.500000 5.331720
25.500000 4.596700
26.500000 3.528202
29.500000 3.058920
27.500000 2.129903
28.500000 1.965610
24.500000 1.563290
e
exit
