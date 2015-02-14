#!/bin/bash
#
# Compile, run the test, plot the output

#------------------------------
prog=test_cardwell

echo "Compiling..."
make $prog

echo "Running program..."
${prog}

#------------------------------
prog=test_cooling

echo "Compiling..."
make $prog

echo "Running program..."
${prog} > ${prog}.out

echo "Plotting..."
gnuplot < plot_${prog}.gp

echo "Displaying"
eog ${prog}.png
