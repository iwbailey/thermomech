#!/usr/bin/gnuplot
# Plot the time vs magnitude from an earthquake catalog

reset

# png
set terminal pngcairo enhanced font 'Verdana,8'
set output 'test_noheat_depthhist.png'

unset key
set grid

set xrange [0:17.5]
set style data histograms
set style fill solid 1.0 border -1

binwidth=17.5/32.0
set boxwidth 0.5*binwidth
bin(x,width)=width*floor(x/width) + binwidth/2.0

set xlabel 'Depth [km]'
set ylabel 'Number Hypocenters'
plot 'test_noheat_eqcat.out' using (bin($3,binwidth)):(1.0) smooth freq with boxes
