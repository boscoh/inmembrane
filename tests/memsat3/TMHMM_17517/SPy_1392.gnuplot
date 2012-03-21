set arrow from 1,1.07 to 6,1.07 nohead lt 3 lw 10
set arrow from 7,1.09 to 29,1.09 nohead lt 1 lw 40
set arrow from 30,1.11 to 43,1.11 nohead lt 4 lw 10
set arrow from 44,1.09 to 66,1.09 nohead lt 1 lw 40
set arrow from 67,1.07 to 72,1.07 nohead lt 3 lw 10
set arrow from 73,1.09 to 93,1.09 nohead lt 1 lw 40
set arrow from 94,1.11 to 102,1.11 nohead lt 4 lw 10
set arrow from 103,1.09 to 125,1.09 nohead lt 1 lw 40
set arrow from 126,1.07 to 136,1.07 nohead lt 3 lw 10
set arrow from 137,1.09 to 159,1.09 nohead lt 1 lw 40
set arrow from 160,1.11 to 163,1.11 nohead lt 4 lw 10
set arrow from 164,1.09 to 186,1.09 nohead lt 1 lw 40
set arrow from 187,1.07 to 220,1.07 nohead lt 3 lw 10
set arrow from 221,1.09 to 243,1.09 nohead lt 1 lw 40
set arrow from 244,1.11 to 257,1.11 nohead lt 4 lw 10
set arrow from 258,1.09 to 280,1.09 nohead lt 1 lw 40
set arrow from 281,1.07 to 286,1.07 nohead lt 3 lw 10
set arrow from 287,1.09 to 306,1.09 nohead lt 1 lw 40
set arrow from 307,1.11 to 310,1.11 nohead lt 4 lw 10
set arrow from 311,1.09 to 330,1.09 nohead lt 1 lw 40
set arrow from 331,1.07 to 342,1.07 nohead lt 3 lw 10
set arrow from 343,1.09 to 365,1.09 nohead lt 1 lw 40
set arrow from 366,1.11 to 374,1.11 nohead lt 4 lw 10
set arrow from 375,1.09 to 394,1.09 nohead lt 1 lw 40
set arrow from 395,1.07 to 398,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for SPy_1392"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:398]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_17517/SPy_1392.eps"
plot "./TMHMM_17517/SPy_1392.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
