#!/usr/bin/gnuplot
# Plot the output from test_loading

reset

ifile = 'test_faultstrength.out'
ofile1 = 'test_faultstrength_1.png'
ofile2 = 'test_faultstrength_2.png'

# png
set terminal pngcairo enhanced font 'Verdana,8'
set output ofile1

# 1:1 ratio
unset key
set size ratio -1

# border
set style line 11 lc rgb '#808080' lt 1
set border 3 front ls 11
set tics nomirror out scale 0.75

# Fault dimensions
set xrange [0:70.0]
set yrange [0:17.5] reverse

set multiplot; # get into multiplot mode
set size 1,0.5; # Each plot uses 1/2 the height

# Plot Current fault strength
set origin 0.0,0.5;
set ylabel 'z [km]'
set cblabel "Strength [MPa]"
set title 'Strength (New)'
plot ifile using 1:2:9 with image

# Plot BZ 1996 result for comparison
set origin 0.0,0.0;
set ylabel 'z [km]'
set xlabel 'x [km]'
set cblabel "Strength [MPa]"
set title 'Strength (BZ1996)'
plot ifile using 1:2:10 with image

unset multiplot # end multiplot mode

print sprintf("written to %s", ofile1)

#-------------------------------------------------------------------------------
# Plot the strain rates
set output ofile2
set multiplot; # get into multiplot mode
set size 1,0.5; # Each plot uses 1/2 the height

# Plot Current fault strength
set origin 0.0,0.5;
set ylabel 'z [km]'
set cblabel "Strain Rate [/yr]"
set title 'Strain Rate at Strength (New)'
plot ifile using 1:2:11 with image

# Plot BZ 1996 result for comparison
set origin 0.0,0.0;
set ylabel 'z [km]'
set xlabel 'x [km]'
set cblabel "Strain Rate [/yr]"
set title 'Strain Rate at Strength (BZ1996)'
plot ifile using 1:2:12 with image

unset multiplot # end multiplot mode

print sprintf("written to %s", ofile2)

# set cbrange [0.001:10]
# set logscale cb
