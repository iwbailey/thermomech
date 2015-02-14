#!/usr/bin/gnuplot
# Plot the output from test_noheat and test_bz1996
reset

# png
set terminal pngcairo enhanced font 'Verdana,8'
set output 'test_noheat_slippartition.png'

unset key
set size ratio -1

# border
set style line 11 lc rgb '#808080' lt 1
set border 3 front ls 11
set tics nomirror out scale 0.75

set xrange [0:128]
set xlabel 'x'

set yrange [0:32] reverse
set ylabel 'z'

set multiplot; # get into multiplot mode
set size 1,0.33;

set origin 0.0,0.66;
set cblabel "Slip From Creep [m]"
plot 'test_noheat_creepslip.txt' using 1:2:3 with image

set origin 0.0,0.33;
set cblabel "Slip From Earthquake [m]"
plot 'test_noheat_eqkslip.txt' using 1:2:3 with image

set origin 0.0,0.0;
set cblabel "Remaining Slip Deficit [m]"
plot 'test_noheat_finalslipdef.txt' using 1:2:3 with image

unset multiplot # exit multiplot mode


# set cbrange [0.001:10]
# set logscale cb
