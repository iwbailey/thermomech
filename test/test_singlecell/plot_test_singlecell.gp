#!/usr/bin/gnuplot
# Plot the output from test_singlecell
reset

# png
set terminal pngcairo enhanced font 'Verdana,8'
set output 'test_singlecell.png'

unset key
#set size ratio -1

# # border
# set style line 11 lc rgb '#808080' lt 1
# set border 3 front ls 11
# set tics nomirror out scale 0.75

#set xrange [0:20]
set xlabel 'time [yr]'


set multiplot; # get into multiplot mode
set size 1,0.33;

# Top plot
set origin 0.0,0.66;
#set yrange [0:32]
set ylabel 'Stress [MPa]'
plot 'test_singlecell.out' using 1:2 w linespoints

set origin 0.0,0.33;
set ylabel 'Temperature [K]'
plot 'test_singlecell.out' using 1:3 w linespoints

set origin 0.0,0.0;
set ylabel 'Creep Velocity [mm/yr]'
plot 'test_singlecell.out' using 1:4 w linespoints

unset multiplot # exit multiplot mode


# set cbrange [0.001:10]
# set logscale cb
