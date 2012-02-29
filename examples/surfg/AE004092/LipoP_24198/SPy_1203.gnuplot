set title "LipoP predictions for SPy_1203"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:5]
set y2range [0:8]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1203.eps"
set arrow from 2,4.32225 to 6,4.32225 nohead lt 1 lw 20
set label "SpI" at 7,4.32225
set arrow from 2,1.70087 to 6,1.70087 nohead lt 4 lw 20
set label "TMH" at 7,1.70087
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,4.32225 to 6,4.32225 nohead lt 1 lw 20
set label "SpI" at 7,4.32225
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
37.500000 7.033730
31.500000 4.460760
33.500000 1.982230
30.500000 0.310670
e
exit
