#!/bin/bash
#
#
# Compile, run the test, plot the output
#
# Requirements:
#  make
#  gnuplot

# What program used to open a png file
if [[ -x /usr/bin/eog ]]; then
    # Look for eog - should be on gnome based systems
    pngviewer=/usr/bin/eog
elif [[ -x /usr/bin/ristretto ]]; then
    # Ristretto should be on xfce based systems
    pngviewer=/usr/bin/ristretto
else
    echo "No png viewer found. Edit the script"
    exit
fi

prog=test_strengthprofile

echo "Compiling..."
make $prog

echo "Running program..."
./${prog} > ${prog}.out

echo "Plotting..."
gnuplot < plot_${prog}.gp

echo "Displaying"
${pngviewer} ${prog}.png


prog=test_faultstrength

echo "Compiling..."
make $prog

echo "Running program..."
./${prog} > ${prog}.out

echo "Plotting..."
gnuplot < plot_${prog}.gp

echo "Displaying"
${pngviewer} ${prog}*.png


prog=test_strength_evol

echo "Compiling..."
make $prog

echo "Running program..."
./${prog} > ${prog}.out

echo "Plotting..."
gnuplot < plot_${prog}.gp

echo "Displaying"
${pngviewer} ${prog}.png
