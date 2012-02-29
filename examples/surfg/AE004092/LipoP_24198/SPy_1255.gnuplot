set title "LipoP predictions for SPy_1255"
set size 2., 1.4
set xrange [1:70] 
set mxtics 10
set yrange [-3:5]
set y2range [0:8]
set ylabel "log-odds"
set term postscript eps color solid "Helvetica" 30
set output "SPy_1255.eps"
set arrow from 2,4.41819 to 6,4.41819 nohead lt 1 lw 20
set label "SpI" at 7,4.41819
set arrow from 2,-0.200913 to 6,-0.200913 nohead lt 3 lw 20
set label "CYT" at 7,-0.200913
set arrow from 2,-2.1601 to 6,-2.1601 nohead lt 4 lw 20
set label "TMH" at 7,-2.1601
set arrow from 2,4.41819 to 6,4.41819 nohead lt 1 lw 20
set label "SpI" at 7,4.41819
# NOTE: The scores below are the log-odds scores with the threshold
# NOTE: subtracted (a hack to make gnuplot make the histogram all
# NOTE: look nice).
plot "-" axes x1y2 title "" with impulses lt 1 lw 20
34.500000 6.271500
36.500000 5.427010
39.500000 5.317290
37.500000 2.751402
38.500000 1.053890
42.500000 0.063120
e
exit
