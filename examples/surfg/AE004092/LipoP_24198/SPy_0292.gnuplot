set title "LipoP predictions for SPy_0292"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:15]
set y2range [0:18]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_0292.eps"
set arrow from 2,11.3779 to 6,11.3779 nohead lt 1 lw 20
set label "SpI" at 7,11.3779
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,11.3779 to 6,11.3779 nohead lt 1 lw 20
set label "SpI" at 7,11.3779
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
22.500000 14.185100
21.500000 11.283430
19.500000 5.936370
20.500000 5.663720
26.500000 5.060340
23.500000 3.751050
16.500000 2.225780
17.500000 1.578660
18.500000 1.160640
24.500000 0.891770
e
exit
