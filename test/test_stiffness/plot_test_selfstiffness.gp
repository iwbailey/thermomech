set term png
set output "test_selfstiffness.png"

set grid back ls 12
set xlabel "Depth [km]"
set ylabel "Stiffness"
set xrange [0:18]
plot './test_selfstiffness.out' using 1:2 with linespoints lw 2 title "Self-stiffness"
