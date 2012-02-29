set title "LipoP predictions for SPy_0337"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0337.eps"
set arrow from 2,6.66808 to 6,6.66808 nohead lt 1 lw 20
set label "SpI" at 7,6.66808
set arrow from 2,2.32928 to 6,2.32928 nohead lt 4 lw 20
set label "TMH" at 7,2.32928
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,6.66808 to 6,6.66808 nohead lt 1 lw 20
set label "SpI" at 7,6.66808
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
40.500000 9.309670
38.500000 7.336000
35.500000 3.395519
37.500000 2.184070
42.500000 0.087690
e
exit
