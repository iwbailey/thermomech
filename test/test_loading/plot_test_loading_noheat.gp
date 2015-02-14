#!/usr/bin/gnuplot
# Plot the output from test_loading

reset

# png
set terminal pngcairo enhanced font 'Verdana,8'
set output 'test_loading_noheat.png'

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

set multiplot; # get into multiplot mode
set size 1,0.5;

set origin 0.0,0.5;
set cblabel "Stress [MPa]"
plot 'test_loading_noheat.out' using 1:2:3 with image

set origin 0.0,0.0;
set cblabel "Slip Deficit [m]"
plot 'test_loading_noheat.out' using 1:2:4 with image
unset multiplot # exit multiplot mode

# set cbrange [0.001:10]
# set logscale cb

reset

# png
set terminal pngcairo enhanced font 'Verdana,8'
set output 'test_loading_noheat_sliprate.png'

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

set cblabel "Creep Rate [mm/yr]"
plot 'test_loading_noheat.out' using 1:2:6 with image
