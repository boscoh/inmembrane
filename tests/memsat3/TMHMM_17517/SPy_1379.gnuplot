set arrow from 1,1.07 to 32,1.07 nohead lt 3 lw 10
set arrow from 33,1.09 to 55,1.09 nohead lt 1 lw 40
set arrow from 56,1.11 to 58,1.11 nohead lt 4 lw 10
set arrow from 59,1.09 to 81,1.09 nohead lt 1 lw 40
set arrow from 82,1.07 to 87,1.07 nohead lt 3 lw 10
set arrow from 88,1.09 to 105,1.09 nohead lt 1 lw 40
set arrow from 106,1.11 to 178,1.11 nohead lt 4 lw 10
set arrow from 179,1.09 to 201,1.09 nohead lt 1 lw 40
set arrow from 202,1.07 to 207,1.07 nohead lt 3 lw 10
set arrow from 208,1.09 to 230,1.09 nohead lt 1 lw 40
set arrow from 231,1.11 to 244,1.11 nohead lt 4 lw 10
set arrow from 245,1.09 to 267,1.09 nohead lt 1 lw 40
set arrow from 268,1.07 to 287,1.07 nohead lt 3 lw 10
set arrow from 288,1.09 to 306,1.09 nohead lt 1 lw 40
set arrow from 307,1.11 to 320,1.11 nohead lt 4 lw 10
set arrow from 321,1.09 to 343,1.09 nohead lt 1 lw 40
set arrow from 344,1.07 to 349,1.07 nohead lt 3 lw 10
set arrow from 350,1.09 to 372,1.09 nohead lt 1 lw 40
set arrow from 373,1.11 to 381,1.11 nohead lt 4 lw 10
set arrow from 382,1.09 to 404,1.09 nohead lt 1 lw 40
set arrow from 405,1.07 to 410,1.07 nohead lt 3 lw 10
set arrow from 411,1.09 to 433,1.09 nohead lt 1 lw 40
set arrow from 434,1.11 to 437,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for SPy_1379"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:437]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_17517/SPy_1379.eps"
plot "./TMHMM_17517/SPy_1379.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
