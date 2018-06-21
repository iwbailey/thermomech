#!/bin/bash
#
# Compile, run the test, plot the output

set -e

prog=test_loading_noheat

echo "Compiling..."
make $prog

echo "Running program..."
./${prog} > ${prog}.out

echo "Plotting..."
gnuplot < plot_${prog}.gp

echo "Displaying"
if [ -x /usr/bin/eog ]; then
    eog ${prog}.png
elif [ -x /usr/bin/ristretto ]; then
    ristretto $prog.png
else
    echo Output check plot written to $prog.png
fi
