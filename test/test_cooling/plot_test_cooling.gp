set term png
set output "test_cooling.png"

set log x
set xrange [0.001:10000]
set grid back ls 12
set xlabel "Time since start of slip [hr]"
set ylabel "Temperature Change [K]"
plot './test_cooling.out' using 1:2 with lines lw 2 title "Slip time = 0.0003 s", \
    './test_cooling.out' using 1:3 with lines lw 2 title "Slip time = 1 hr", \
    './test_cooling.out' using 1:4 with lines lw 2 title "Slip time = 1 day"

