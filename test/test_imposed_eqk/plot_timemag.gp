#!/usr/bin/gnuplot
# Plot the time vs magnitude from an earthquake catalog

reset

# png
set terminal pngcairo enhanced font 'Verdana,8'
set output 'test_circlerupt_timemag.png'

unset key
set grid
set size ratio 0.5

set xrange [0:1.0]
set xlabel 'Time [yr]'
set ylabel 'M_W'


plot 'eqcat.out' using 1:5 with points pt 7 ps 0.5
