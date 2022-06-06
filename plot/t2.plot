set terminal postscript eps color
set output 't2.eps'
unset key


plot 'explicit.dat' u 1:3 with lines lw 3 title 't=0'

