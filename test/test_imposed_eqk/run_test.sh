#!/bin/bash
#
# Compile, run the test, plot the output

#------------------------------
prog=test_circlerupt

echo "Compiling..."
make $prog

echo "Running program..."
./${prog} > eqcat.out


for plotname in magfreq hypoxy timemag slippartition stresschange tempchange; do
    echo "Plotting ${plotname}..."
    gnuplot < plot_${plotname}.gp
    echo "Written to ${prog}_${plotname}.png"
done
