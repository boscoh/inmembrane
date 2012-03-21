set arrow from 1,1.11 to 14,1.11 nohead lt 4 lw 10
set arrow from 15,1.09 to 30,1.09 nohead lt 1 lw 40
set arrow from 31,1.07 to 36,1.07 nohead lt 3 lw 10
set arrow from 37,1.09 to 59,1.09 nohead lt 1 lw 40
set arrow from 60,1.11 to 91,1.11 nohead lt 4 lw 10
set arrow from 92,1.09 to 114,1.09 nohead lt 1 lw 40
set arrow from 115,1.07 to 126,1.07 nohead lt 3 lw 10
set arrow from 127,1.09 to 149,1.09 nohead lt 1 lw 40
set arrow from 150,1.11 to 179,1.11 nohead lt 4 lw 10
set arrow from 180,1.09 to 197,1.09 nohead lt 1 lw 40
set arrow from 198,1.07 to 216,1.07 nohead lt 3 lw 10
set arrow from 217,1.09 to 239,1.09 nohead lt 1 lw 40
set arrow from 240,1.11 to 258,1.11 nohead lt 4 lw 10
set arrow from 259,1.09 to 281,1.09 nohead lt 1 lw 40
set arrow from 282,1.07 to 301,1.07 nohead lt 3 lw 10
set arrow from 302,1.09 to 324,1.09 nohead lt 1 lw 40
set arrow from 325,1.11 to 333,1.11 nohead lt 4 lw 10
set arrow from 334,1.09 to 356,1.09 nohead lt 1 lw 40
set arrow from 357,1.07 to 376,1.07 nohead lt 3 lw 10
set arrow from 377,1.09 to 399,1.09 nohead lt 1 lw 40
set arrow from 400,1.11 to 411,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for SPy_1949"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:411]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_17517/SPy_1949.eps"
plot "./TMHMM_17517/SPy_1949.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
