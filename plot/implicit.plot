set terminal postscript eps color
set output 'implicit.eps'
set xrange[0:1]
set xtics 0,0.2,1.0 font ",20"
set mxtics 5
set xlabe "x" font ",25" offset 0,1,0
set yrange[0.:1]
set ytics 0,0.2,1 font ",20"
set mytics 5
set ylabel "u" font ",25" offset 1,-0.5,0
set key nobox
set key spacing 2.0 
#set key at 0.9,4.5
set border lw 1.5
set bmargin at screen 0.1
set tmargin at screen 0.92
#set title 'implicit' font ',30' offset 0,-1,0
plot 'implicit.dat' u 1:2 with lines lw 3 title 'initial','implicit.dat' u 1:3 pt 7 ps 0.5 title 'numerical solution','implicit.dat' u 1:4 with lines lw 3 title 'exact solution' 