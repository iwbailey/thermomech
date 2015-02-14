#!/usr/bin/gnuplot
# Plot the magnitude frequency distribution from an earthquake catalog

reset

# png
set terminal pngcairo enhanced font 'Verdana,8'
set output 'test_noheat_magfreq.png'

unset key
set grid
set size ratio 1

set xrange [3:7]
set xlabel 'M_W Exceeded'
set ylabel 'Frequency [/yr]'
set logscale y

cumm_sum=0.0
nrec=0.0
CDF(x)=(nrec=nrec+1.0/100, nrec)

#plot 'test_bz1996.out' using 5:(-1.0) smooth cumul
plot '< sort -k5r test_noheat_eqcat.out' using 5:(CDF($1)) w lines linewidth 2
