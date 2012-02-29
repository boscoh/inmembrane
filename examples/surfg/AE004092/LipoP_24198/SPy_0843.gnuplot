set title "LipoP predictions for SPy_0843"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0843.eps"
set arrow from 2,8.07874 to 6,8.07874 nohead lt 1 lw 20
set label "SpI" at 7,8.07874
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,8.07874 to 6,8.07874 nohead lt 1 lw 20
set label "SpI" at 7,8.07874
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
21.500000 10.988750
26.500000 5.975960
19.500000 5.333650
20.500000 3.314870
24.500000 3.045128
22.500000 2.532802
17.500000 0.263440
18.500000 0.113090
e
exit
