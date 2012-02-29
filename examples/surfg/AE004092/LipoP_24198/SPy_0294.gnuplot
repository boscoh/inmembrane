set title "LipoP predictions for SPy_0294"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:5]
set y2range [0:8]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0294.eps"
set arrow from 2,3.72404 to 6,3.72404 nohead lt 1 lw 20
set label "SpI" at 7,3.72404
set arrow from 2,3.06046 to 6,3.06046 nohead lt 4 lw 20
set label "TMH" at 7,3.06046
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,3.72404 to 6,3.72404 nohead lt 1 lw 20
set label "SpI" at 7,3.72404
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
29.500000 6.157620
27.500000 3.869622
31.500000 2.239926
25.500000 2.166062
23.500000 2.135585
30.500000 1.389300
21.500000 0.326110
e
exit
