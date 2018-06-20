#!/bin/bash
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

# Run the stiffness test
prog=test_stiffness

echo "Compiling ${prog}..."
make $prog

echo "Running program ${prog}..."
./${prog} > ${prog}.out

echo "Plotting ${prog}..."
gnuplot < plot_${prog}.gp

echo "Displaying ${prog}..."
${pngviewer} ${prog}.png

# Run the self stiffness test
prog=test_selfstiffness

echo "Compiling ${prog}..."
make $prog

echo "Running program ${prog}..."
./${prog} > ${prog}.out

echo "Plotting ${prog}..."
gnuplot < plot_${prog}.gp

echo "Displaying ${prog}..."
${pngviewer} ${prog}.png
