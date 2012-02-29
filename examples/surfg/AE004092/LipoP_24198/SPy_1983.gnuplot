set title "LipoP predictions for SPy_1983"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:15]
set y2range [0:18]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1983.eps"
set arrow from 2,11.3534 to 6,11.3534 nohead lt 1 lw 20
set label "SpI" at 7,11.3534
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,11.3534 to 6,11.3534 nohead lt 1 lw 20
set label "SpI" at 7,11.3534
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
37.500000 14.179800
36.500000 10.797590
34.500000 7.615390
35.500000 7.606860
32.500000 7.399180
33.500000 3.867878
39.500000 3.434973
29.500000 2.207482
31.500000 0.568930
27.500000 0.205290
e
exit
