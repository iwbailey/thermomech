#!/bin/bash
#
# Compile, run the test, plot the output

# Run the stiffness test
prog=test_stiffness

echo "Compiling ${prog}..."
make $prog

echo "Running program ${prog}..."
${prog} > ${prog}.out

echo "Plotting ${prog}..."
gnuplot < plot_${prog}.gp

echo "Displaying ${prog}..."
eog ${prog}.png

prog=test_selfstiffness

echo "Compiling ${prog}..."
make $prog

echo "Running program ${prog}..."
${prog} > ${prog}.out

echo "Plotting ${prog}..."
gnuplot < plot_${prog}.gp

echo "Displaying ${prog}..."
eog ${prog}.png

