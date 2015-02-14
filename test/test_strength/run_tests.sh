#!/bin/bash
#
# Compile, run the test, plot the output

prog=test_strengthprofile

echo "Compiling..."
make $prog

echo "Running program..."
${prog} > ${prog}.out

echo "Plotting..."
gnuplot < plot_${prog}.gp

echo "Displaying"
eog ${prog}.png


prog=test_faultstrength

echo "Compiling..."
make $prog

echo "Running program..."
${prog} > ${prog}.out

echo "Plotting..."
gnuplot < plot_${prog}.gp

echo "Displaying"
eog ${prog}*.png


prog=test_strength_evol

echo "Compiling..."
make $prog

echo "Running program..."
${prog} > ${prog}.out

echo "Plotting..."
gnuplot < plot_${prog}.gp

echo "Displaying"
eog ${prog}.png
