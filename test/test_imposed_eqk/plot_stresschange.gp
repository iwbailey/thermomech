#!/usr/bin/gnuplot
# Plot the output from test_loading

reset

# png
set terminal pngcairo enhanced font 'Verdana,8'
set output 'test_circlerupt_stresschange.png'

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
set size 1,0.35;

set origin 0.0,0.6;
set cblabel "Initial Stress [MPa]"
plot 'test_circlerupt_initstress.txt' using 1:2:3 with image

set origin 0.0,0.3;
set cblabel "Final Stress [MPa]"
plot 'test_circlerupt_finalstress.txt' using 1:2:3 with image

set origin 0.0,0.0;
set cblabel "Stress Increase [MPa]"
plot '< paste test_circlerupt_initstress.txt test_circlerupt_finalstress.txt' \
     using 1:2:($6-$3) with image

unset multiplot # exit multiplot mode

# set cbrange [0.001:10]
# set logscale cb
