#!/usr/bin/gnuplot
# Plot the hypocenter locations on the fault

reset

# png
set terminal pngcairo enhanced font 'Verdana,8'
set output 'test_circlerupt_hypoxy.png'

unset key
set grid
set size ratio -2

set xrange [0:70.0]
set yrange [0:17.5] reverse
set xlabel 'x [km]'
set ylabel 'y [km]'
set cblabel "Time [yr]"

plot 'eqcat.out' using 2:3:(1.25*$5-3.5):1 with points pt 6 ps variable \
     lt palette
