#!/usr/bin/gnuplot
# Plot the output from test_stiffness

reset

# # wxt
# set terminal wxt size 350,262 enhanced font 'Verdana,10' persist

# png
set terminal pngcairo enhanced font 'Verdana,8'
set output 'test_stiffness.png'

unset key
set size ratio -1

# border
set style line 11 lc rgb '#808080' lt 1
set border 3 front ls 11
set tics nomirror out scale 0.75

set xrange [0:70.0]
set xlabel 'x [km]'

set yrange [0:17.5] reverse
set ylabel 'z [km]'

set cbrange [0.001:10]
set logscale cb
set cblabel "Stress change [Mpa]


plot 'test_stiffness.out' using 1:2:3 with image
