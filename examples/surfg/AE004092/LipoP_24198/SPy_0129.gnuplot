set title "LipoP predictions for SPy_0129"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:10]
set y2range [0:13]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0129.eps"
set arrow from 2,5.57386 to 6,5.57386 nohead lt 4 lw 20
set label "TMH" at 7,5.57386
set arrow from 2,1.21971 to 6,1.21971 nohead lt 1 lw 20
set label "SpI" at 7,1.21971
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,5.57386 to 6,5.57386 nohead lt 4 lw 20
set label "TMH" at 7,5.57386
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
34.500000 2.267271
29.500000 2.050258
35.500000 1.805350
33.500000 0.732320
28.500000 0.340880
e
exit
