set title "LipoP predictions for SPy_2010"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:15]
set y2range [0:18]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_2010.eps"
set arrow from 2,12.0873 to 6,12.0873 nohead lt 1 lw 20
set label "SpI" at 7,12.0873
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,12.0873 to 6,12.0873 nohead lt 1 lw 20
set label "SpI" at 7,12.0873
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
27.500000 14.196200
25.500000 13.886800
30.500000 8.483470
29.500000 8.301890
24.500000 6.592050
23.500000 6.041620
26.500000 4.987080
31.500000 4.153320
22.500000 2.438279
33.500000 0.075430
e
exit
