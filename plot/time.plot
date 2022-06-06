set terminal postscript eps color
set output 'time.eps'

set xrange[0:8]
set xtics 0,2,8 font ",20"
set mxtics 5
set xlabe "N" font ",25" offset 0,0,0

set mytics 5
set ylabe "t" font ",25" offset 0,0,0

set key box lw 2.0 font ',20'
set key spacing 1.5 
set key left
set border lw 1.5
set bmargin at screen 0.1
set tmargin at screen 0.92
#set title 'Time cost' font ',30' offset 0,-1,0

plot 'time.dat' u 1:2 with lp lw 3 title 'explicit', 'time.dat' u 1:3 with lp lw 3 title 'implicit'