set terminal postscript eps color
set output 'flops.eps'

set xrange[0:8]
set xtics 0,2,8 font ",20"
set mxtics 5
set xlabe "N" font ",25" offset 0,0,0

set yrange[1.0e4:2.0e7]
set ytics  font ",20"
set mytics 5
set ylabe "flops/sec" font ",25" offset 0,0,0

set key box lw 2.0 font ',20'
set key spacing 1.5 
set key right
set border lw 1.5
set bmargin at screen 0.1
set tmargin at screen 0.92

plot 'flops.dat' u 1:2 with lp lw 3 title 'explicit', 'flops.dat' u 1:3 with lp lw 3 title 'implicit'