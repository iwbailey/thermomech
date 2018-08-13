#!/bin/bash
#
# Compile, run the test, plot the output

#------------------------------
prog=test_cardwell

echo "Compiling..."
make $prog

echo "Running program..."
./${prog}

#------------------------------
prog=test_cooling

echo "Compiling..."
make $prog

echo -e "\n\nRunning program..."
./${prog} > ${prog}.out

echo "Plotting..."
gnuplot < plot_${prog}.gp

echo "Written output to" ${prog}.png
