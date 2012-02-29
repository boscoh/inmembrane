set title "LipoP predictions for SPy_1357"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1357.eps"
set arrow from 2,4.8635 to 6,4.8635 nohead lt 1 lw 20
set label "SpI" at 7,4.8635
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,4.8635 to 6,4.8635 nohead lt 1 lw 20
set label "SpI" at 7,4.8635
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
33.500000 7.726860
32.500000 3.475991
34.500000 2.641895
35.500000 0.391430
27.500000 0.359530
e
exit
