set title "LipoP predictions for SPy_0031"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:5]
set y2range [0:8]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0031.eps"
set arrow from 2,3.21838 to 6,3.21838 nohead lt 1 lw 20
set label "SpI" at 7,3.21838
set arrow from 2,1.9338 to 6,1.9338 nohead lt 4 lw 20
set label "TMH" at 7,1.9338
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,3.21838 to 6,3.21838 nohead lt 1 lw 20
set label "SpI" at 7,3.21838
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
25.500000 4.675230
30.500000 4.338120
23.500000 3.933658
22.500000 2.105246
21.500000 1.961510
29.500000 0.803730
28.500000 0.717510
e
exit
